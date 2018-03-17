%% Multi-scale Dictionary Learning (simulation)
%
% This m.file operates its calculation based on
% learning_recurrent_waveforms_v1. You can get the scripts from reference.
%
% reference: A. J. Brockmeier et al., 2016, "Learning Recurrent Waveforms within EEGs"
%
% Nested folder and function:
% - decompose
%   - blockbasedMP
%   - multipassMP
%   - multistaegMP
%   - shiftinvariantMP
% - learning
%   - MPSVD
%   - multistageWaveformLearning
% - utilities
%   - atomsToSources
%   - fasterXcorr
%   - sourcesToSignalComponents
%
% Suguru Kanoga, last modification 27 Sep. 2017
% s.kanouga@aist.go.jp

function [results] = sep_MSDL(input_data,correct_eeg,testing_artifact,Fs,time_length,code_dir)

addpath decompose/
addpath learning/
addpath utilities/

% parameters for multi-scale dictionary learning
coverage = 1;
learning_approximation_passes = 4;
do_nonneg = 1;

coeff_per_filter = [16 32 64 125 250];
block_start_indices = [];
signal_params = struct('block_starts',block_start_indices,'freq',Fs);
subsample_rates = ones(1,size(coeff_per_filter,2)); % decimation rates
nfilt = 4;
Filters_per_scale = nfilt*ones(size(coeff_per_filter));
shared_model_params = {'model_sizes',cat(2,coeff_per_filter(:),Filters_per_scale(:)),...
                       'subsample_rates',subsample_rates,'freq',Fs,...
                       'coverage',coverage,...
                       'nonneg_flag',do_nonneg,...
                       'approximation_passess',learning_approximation_passes};
N = size(input_data,2);

d = designfilt('bandpassiir',...
               'FilterOrder',4,...
               'HalfPowerFrequency1',4,...
               'HalfPowerFrequency2',40,...
               'DesignMethod','butter',...
               'SampleRate',Fs);

for sn_variation = 1:size(input_data,3)
    re_eeg_ep = squeeze(input_data(:,:,sn_variation));
    re_artifact_ep = zeros(size(re_eeg_ep));
    for ep = 1:N
        x = squeeze(input_data(:,ep,sn_variation));
        amp_thr = mean(x)+2*std(x);
        freq_thresh = 3/numel(x);
        c = cat(2,shared_model_params,{'minimum_freq',freq_thresh});
        model_params = cell2struct(c(2:2:end),c(1:2:end),2);

        % learning new dictionaries incluging ocular artifactual dictionaries
        tic
        [MSDict,~,~] = multistageWaveformLearning(x,signal_params,model_params);
        results.MSDict(ep,sn_variation) = {MSDict};

        group_dim = cellfun(@(x) size(x,2), MSDict);

        atoms_all = multistageMP(x,MSDict,signal_params,model_params);
        Natoms = sum(cellfun(@(x) size(x,1), atoms_all));
        estimate_sources = atomsToSources(atoms_all,group_dim,numel(x),Natoms);
        dammy_estimate_sources = estimate_sources;
        
        for ii = 1:length(coeff_per_filter)
            for jj = 1:size(estimate_sources{1,ii},2)
                [I,~,W] = find(estimate_sources{1,ii}(:,jj));
                ind_lib = zeros(1,size(W,1)); oa_ind_lib = zeros(1,size(W,1));
                for kk = 1:size(W,1)
                    temp_source = estimate_sources{1,ii}(:,jj);
                    dammy_temp_source = dammy_estimate_sources{1,ii}(:,jj);
                    % identification
                    if sum(ge(abs(full(temp_source(I(kk,:))*MSDict{ii,:}(:,jj))),amp_thr)) >= 1
                       ind_lib(:,kk) = kk; 
                    else
                        oa_ind_lib(:,kk) = kk;
                    end
                end
                ind_lib = find(ind_lib);
                oa_ind_lib = find(oa_ind_lib);
                
                dammy_temp_source(I(oa_ind_lib,:)) = zeros;
                temp_source(I(ind_lib,:)) = zeros;
                
                estimate_sources{1,ii}(:,jj) = temp_source;
                dammy_estimate_sources{1,ii}(:,jj) = dammy_temp_source;
            end
        end
        
        [VV_eeg,~] = sourcesToSignalComponents(MSDict,estimate_sources);
        re_eeg_ep(:,ep) = full(sum(cell2mat(VV_eeg(:)'),2));           
        re_eeg_ep(:,ep) = filtfilt(d,re_eeg_ep(:,ep));
        
        results.ProcessingTime(ep,sn_variation) = {toc};     
                
        [VV_oa,~] = sourcesToSignalComponents(MSDict,dammy_estimate_sources);
        re_artifact_ep(:,ep) = full(sum(cell2mat(VV_oa(:)'),2)); 
        re_artifact_ep(:,ep) = filtfilt(d,re_artifact_ep(:,ep));
        results.MSDict(ep,sn_variation) = {MSDict};  
                
        %%%%%%%%%%%%%%
        % evaluation %
        %%%%%%%%%%%%%%
        results.mae_eeg(ep,sn_variation) = {sum(abs(correct_eeg(:,ep)-re_eeg_ep(:,ep)))/(Fs*time_length)};
        results.mae_artifact(ep,sn_variation) = {sum(abs(testing_artifact(:,ep)-re_artifact_ep(:,ep)))/(Fs*time_length)};
    end
    
    results.re_eeg(sn_variation) = {re_eeg_ep};
    results.re_artifact(sn_variation) = {re_artifact_ep};
end

cd(code_dir);
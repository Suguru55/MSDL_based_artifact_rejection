%% Semi-simulated data making (main)
%
% This m.file makes semi-simulated data from dataset IIa of BCI
% competition IV (http://www.bbci.de/competition/iv/#datasets).
% 
% [1] Epcohs are segmented through trigg function
% [2] Artifact-contaminated epochs (HDR.ArtifactSelection = 1) are extracted and decomposed into
%     ICs by InfoMax ICA implemented in EEGLAB (https://sccn.ucsd.edu/eeglab/download.php)
% [3] The decomposed ICs are identified as EEG or ocular artifact (OA) by
%     soft thresholds based on indecies of modified sample entropy and kurtosis
% [4] The identified artifactual ICs are further decomposed into subbands
%     through discrete wavelet transform and then an additional soft threshold performs to identify EEG or artifactual subband.
% [5] The identified artifactual subbands are going to reconstruct by inverse wavelet transform
% [6] Inverse wavelet transform and linear processing of ICA will extract independent OA sources S_OA from artifactual epochs
% [7] The sources are superimposed on the C3 EEG data (HDR.ArtifactSelection = 0) as follows:
%        X = X_EEG + lambda * S_OA,
%     with lambda representing the contribution of OA activity.
%     The superimposed EEG epochs will be randomly chosen from all EEG
%     epochs (we make 20 validation combinations).
%     SNR can be adjusted by changing the parameter lambda:
%        SNR = RMS(X_EEG)/RMS(lambda*S_OA),
%     where RMS value is defined as:
%        RMS(X) = sqrt(sum(X(c,t)^2)/(C*T)),
%     where C and T are the number of EEG channel and the number of
%     samples.
% [8] SNR values from 0.5 to 1.5 with a 0.1 step are prepared. 
%
%  !! We use only class 1 and 2 for classification (MI_classification.m)
%
%  !! Each subject has different number of artifactual epochs (A06T is having many artifactual epochs)
%
% Nested functions:
% - automatic_wICA
%   - runica
%
% Needed toolboxes:
% - Statics and Machine Learning Toolbox
% - Wavelet Toolbox
% 
% Suguru Kanoga, last modification 19 Oct. 2017
% s.kanouga@aist.go.jp

clear all
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% config
Fs = 250;         % sampling rate
time_length = 6;  % 6-s epochs
sub_num = 9;
ch_all = 25;
ch_eeg = 22;
fh = 40;          % band-pass frequency range
fl = 4;
position = 8;     % C3
var_size = 20;    % 20 variation datasets
sn = 0.5:0.1:1.5; % SNR values from 0.5 to 1.5 with a 0.1 step

thre_order = [1,3,6; 3,11,12; 3,4,7; 14,15,22; 1,7,9;...
    2,11,36; 2,9,17; 13,16,19; 6,8,14];    % order index for data_driven threshold

% change this directories to yours
data_dir = 'C:\Users\sdxwh\Desktop\MSDL_based_artifact_rejection-master\data';
code_dir = 'C:\Users\sdxwh\Desktop\MSDL_based_artifact_rejection-master';
addpath('C:\toolbox\biosig4octmat-3.1.0\biosig\t250_ArtifactPreProcessingQualityControl')    % to use trigg function
addpath('C:\toolbox\biosig4octmat-3.1.0\biosig\t200_FileAccess');                            % to use sload function
% !! I guess BioSig toolbox has specialized function for solving optimal
%    problem that deteriorates performance of lassoglm, so I avoided
%    storing full paths.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract artifact-contaminated epochs from training dataset
for sub_id = 1:sub_num
    %%%%%%%%%%%%%
    % load data %
    %%%%%%%%%%%%%
    cd(data_dir);
        
    eval(sprintf('filename = [''A0%dT.gdf'']',sub_id));
    [s, HDR] = sload(filename, 0);
    eval(sprintf('filename = [''A0%dT_label.mat''];',sub_id));
    load(filename);
        
    %%%%%%%%%%%%
    % fill NaN %
    %%%%%%%%%%%%
    s = fillmissing(s,'spline');
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % filtering (4-40 Hz) %
    %%%%%%%%%%%%%%%%%%%%%%%
    % 4th-order butterworth filter
    d = designfilt('bandpassiir','FilterOrder',4,...
                   'HalfPowerFrequency1',fl,'HalfPowerFrequency2',fh,...
                   'DesignMethod','butter','SampleRate',Fs);
    
    f_s = zeros(size(s));
    
    for i = 1:size(s,2)
        f_s(:,i) = filtfilt(d,s(:,i)); 
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sprit continuous signal into epochs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % you should install BioSig toolbox (http://biosig.sf.net/)
    % but don't add the path constantly because this toolbox changes
    % results of L1 norm regularization (I don't know the detail...)
    % that decreases the performance of lasso logistic regression
    PRE = 0;
    PST = time_length*Fs-1;
    
    [X,sz] = trigg(f_s, HDR.TRIG, PRE, PST);
    epochs = reshape(X,sz);
    cd(code_dir);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract artifactual epochs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind = find(HDR.ArtifactSelection); % provided by BCI competition
    
    artifact_epochs = epochs(:,:,ind);
    dammy_classlabel = classlabel;
    
    epochs(:,:,ind) = [];
    dammy_classlabel(ind,:) = [];
    
    eeg_epochs = epochs;
    eeg_label = dammy_classlabel;    
    
    % !! change the order of artifactual epochs
    %    Å® this process will reduce computational cost for source
    %       identification and improve approximation quality in ICA
    
    re_ind = thre_order(sub_id,:);
    re_order = 1:size(ind,1);
    re_order = re_order';
    re_order(re_ind,:) = [];
    re_order = [re_ind';re_order];
    
    artifact_sig = reshape(artifact_epochs(:,:,re_order),[sz(1,1),sz(1,2)*size(ind,1)]);
    artifact_id = ind(re_order,:);
    artifact_label = classlabel(ind(re_order,:),:);

    cd(data_dir);
    eval(sprintf('filename = [''A0%d_Epochs.mat'']',sub_id));
    save(filename,'artifact_sig','artifact_label','artifact_id','sz',...
                  'eeg_epochs','eeg_label');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal separation based on wICA (Infomax ICA algorithm)
time = 18;   % time window to calculate modified multiscale sample entropy
             % first 3 epochs (4500 data points) are used for this identification
             %    When you use the index of sample entropy, the data points should
             %    be ranged from 100 to 5000.
             %    references: [1] D.Lake et al., 2002, "Sample Entropy Analysis of Neonatal Heart Rate Varibility"
             %                [2] S. M. Pincus, 1991, "Approximate Entropy as a Measure of System Complexity"
             %                [3] D. Wand et al., 2012, "Multi-class Motor Imagery EEG Decoding for Brain-computer Intefaces"

for sub_id = 1:sub_num
    %%%%%%%%%%%%%
    % load data %
    %%%%%%%%%%%%%
    cd(data_dir);
    
    eval(sprintf('filename = [''A0%d_Epochs.mat'']',sub_id));
    load(filename);
    input_data = artifact_sig;

    cd(code_dir);
    [artifact,eeg,properties.weight,properties.icasig,properties.ICnum] = automatic_wICA(input_data,Fs,ch_all,time);
    
    eeg_res = eeg;
    artifact_res = artifact;
    
    input_data = reshape(input_data,[sz(1,1:2) size(artifact_id,1)]);
    artifact = reshape(artifact,[sz(1,1:2) size(artifact_id,1)]);
    eeg = reshape(eeg,[sz(1,1:2) size(artifact_id,1)]);
        
    cd(data_dir);
    eval(sprintf('filename = [''A0%d_Separated.mat'']',sub_id));
    save(filename,'artifact','eeg','properties','input_data');  
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split artifactual signal again based on class labels for testing and
% training data of artifact rejection methods

for sub_id = 1:sub_num
    cd(data_dir);
    eval(sprintf('filename = [''A0%d_Epochs.mat'']',sub_id));
    load(filename);
    eval(sprintf('filename = [''A0%d_Separated.mat'']',sub_id));
    load(filename);
    
    eeg_epochs_1 = eeg_epochs(:,:,eeg_label==1);
    eeg_epochs_2 = eeg_epochs(:,:,eeg_label==2);
    
    testing_artifact = artifact(position(1,1),:,artifact_label<3); % X_OA
    sim_artifact_label = artifact_label(artifact_label<3);
    
    eog_channel = input_data(23:25,:,artifact_label<3);     % reference channels for adaptive filter
    testing_artifact = cat(1,testing_artifact,eog_channel);
    
    rms_oa = sqrt(sum(squeeze(testing_artifact(1,:,:)).^2)./(Fs*time_length));
    
    training_artifact = artifact(:,:,artifact_label>2); % for supervised algotirhms
    
    var_ind = zeros(var_size,sum(artifact_label<3));
    for iter = 1:var_size
        var_1 = randperm(size(eeg_epochs_1,3),sum(artifact_label==1));
        var_2 = randperm(size(eeg_epochs_2,3),sum(artifact_label==2));
        a1 = 1; a2 = 1;
        target_eeg = zeros(Fs*time_length,size(sim_artifact_label,1));
        for i = 1:size(sim_artifact_label,1)
            % X_EEG
            if sim_artifact_label(i,:) == 1
                target_eeg(:,i) = squeeze(eeg_epochs_1(position(1,1),:,var_1(:,a1))); % class 1 X_EEG
                var_ind(iter,i) = var_1(:,a1);
                a1 = a1 + 1;
            else
                target_eeg(:,i) = squeeze(eeg_epochs_2(position(1,1),:,var_2(:,a2))); % class 2 X_EEG
                var_ind(iter,i) = var_2(:,a2); 
                a2 = a2 + 1;               
            end
        end
        rms_eeg = sqrt(sum(target_eeg.^2)./(Fs*time_length));
        lambda = zeros(size(sn,2),size(target_eeg,2));
        X = zeros(Fs*time_length,size(target_eeg,2),size(sn,2));

        for sn_variation = 1:size(sn,2)
            SN = sn(:,sn_variation);
            lambda(sn_variation,:) = rms_eeg./(rms_oa*SN);
            X(:,:,sn_variation) = target_eeg + lambda(sn_variation,:).*squeeze(testing_artifact(1,:,:));
        end
        
        if iter == 1
            X_sets = X;
            targets = target_eeg;
            lambdas = lambda;
        else
            X_sets = cat(2,X_sets,X);
            targets = cat(2,targets,target_eeg);
            lambdas = cat(2,lambdas,lambda);
        end
    end

    cd(data_dir);
    eval(sprintf('filename = [''A0%d_SS_data.mat'']',sub_id));
    save(filename,'var_ind','testing_artifact','training_artifact','targets','X_sets','lambdas');  
end
cd(code_dir);

% check
% figure(4)
% for i = 1:size(target_eeg,2)
% plot(target_eeg(:,i),'g--','LineWidth',1.5);hold on
% plot(squeeze(testing_artifact(1,:,i)),'r--','LineWidth',1.5);
% plot(X(:,i,1),'k','LineWidth',2);
% plot(X(:,i,11),'k','LineWidth',1.5);hold off
% grid on
% set(gca,'FontName','Times New Roman','FontSize',20)
% xlabel('Sample')
% ylabel('Amplitude (\muV)')
% ylim([-100 100])
% legend('Clean EEG','Ocular artifact','Noisy EEG (SNR=0.5)','Noisy EEG (SNR=1.5)');
% pause
% end
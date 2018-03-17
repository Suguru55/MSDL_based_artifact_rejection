%% Single channel artifact reductionn (main) 
%
% This m.file handles artifact reduction and evaluates waveform
% representation ability using semi-simulated data
%
% Before the use of this code, you should make semi-simulation data by
% Make_ss_data.m
%
% Suguru Kanoga, last modification 2 Oct. 2017
% s.kanouga@aist.go.jp

clear all
close all

% change your directory
data_dir = 'C:\Users\sdxwh\Desktop\MSDL\data';
code_dir = 'C:\Users\sdxwh\Desktop\MSDL';
save_dir = 'C:\Users\sdxwh\Desktop\MSDL';
addpath(genpath('C:\Users\sdxwh\Desktop\MSDL\minFunc_2012')) % add all sub-directories to the path
addpath(genpath('C:\Users\sdxwh\Desktop\MSDL\MIToolbox-3.0.0\matlab'))

% config
Fs = 250;
sub_num = 9;
position = 8;    % C3
time_length = 6; % 6-s epoch for evaluation
var_size = 20;   % 20 variation datasets
ref_ch = [2,3];  % reference channels for adaptive filter (ch #1 and #2)

% This 'for' loop takes over 1 day, so I recommend selecting one subject and
% one iteration for checking the performance (Even if you selected that setting, it may take about 10 minutes...).
% ex.) sub_id = 1, iter = 1
for sub_id = 1:sub_num
    % load data  
    cd(data_dir);
    eval(sprintf('filename = [''A0%d_Epochs.mat'']',sub_id));
    load(filename);
    eval(sprintf('filename = [''A0%d_SS_Data.mat'']',sub_id));
    load(filename);
    eval(sprintf('filename = [''A0%d_Separated.mat'']',sub_id));
    load(filename);

    for iter = 1:var_size
        clear -global lam_num
        input_data = X_sets(:,1+(iter-1)*size(var_ind,2):size(var_ind,2)+(iter-1)*size(var_ind,2),:);
        correct_eeg = targets(:,1+(iter-1)*size(var_ind,2):size(var_ind,2)+(iter-1)*size(var_ind,2));
    
        % artifact reduction (signal separation)
        cd(code_dir);
        AF(iter) = sep_AF(input_data,correct_eeg,squeeze(testing_artifact(1,:,:)),testing_artifact(ref_ch,:,:),Fs,time_length,code_dir);
        DOAC(iter) = sep_DOAC(input_data,correct_eeg,squeeze(testing_artifact(1,:,:)),Fs,time_length,code_dir,save_dir,artifact_label(artifact_label<3,:)); 
        % DOAC needs many time for training
        EMDandCC(iter) = sep_EMDandCC(input_data,correct_eeg,squeeze(testing_artifact(1,:,:)),Fs,time_length,code_dir);
        MSDL(iter) = sep_MSDL(input_data,correct_eeg,squeeze(testing_artifact(1,:,:)),Fs,time_length,code_dir);
        fprintf('Sub %i iteration %i\n',sub_id,iter);
        % show MAE for artifact-reduced EEG data when SNR=0.5
        fprintf('AF: %4.2f, DOAC: %4.2f, EMD-CC: %4.2f, MSDL: %4.2f',...
            mean(cell2mat(AF.mae_eeg(:,1))),mean(cell2mat(DOAC.mae_eeg(:,1))),mean(cell2mat(EMDandCC.mae_eeg(:,1))),mean(cell2mat(MSDL.mae_eeg(:,1))));
    end 
end

% Visualization(one trial)
figure(1)
subplot(2,1,1)
plot(AF(1).re_eeg{1,1}(:,1),'y','LineWidth',2);hold on
plot(DOAC(1).re_eeg{1,1}(:,1),'g','LineWidth',2);
plot(EMDandCC(1).re_eeg{1,1}(:,1),'r','LineWidth',2);
plot(MSDL(1).re_eeg{1,1}(:,1),'b','LineWidth',2);
plot(correct_eeg(:,1),'m--');
plot(input_data(:,1,1),'k-.');hold off
grid on
set(gca,'FontName','Times new roman','FontSize',16)
legend('X\ast_{EEG} (AF)','X\ast_{EEG} (DOAC)','X\ast_{EEG} (EMD-CC)','X_{EEG}(MSDL)','X_{EEG}','X_{SIM} (SNR = 0.5)')
ylabel('Amplitude (\muV)')
ylim([-70 70])
xlim([1 Fs*time_length])
subplot(2,1,2)
plot(AF(1).re_artifact{1,1}(:,1),'y','LineWidth',2);hold on
plot(DOAC(1).re_artifact{1,1}(:,1),'g','LineWidth',2);
plot(EMDandCC(1).re_artifact{1,1}(:,1),'r','LineWidth',2);
plot(MSDL(1).re_artifact{1,1}(:,1),'b','LineWidth',2);    
plot(input_data(:,1,1)-correct_eeg(:,1),'m--');
plot(input_data(:,1,1),'k-.');hold off
grid on
set(gca,'Fontname','Times New Roman','Fontsize',16)
legend('X\ast_{OA} (AF)','X\ast_{OA} (DOAC)','X\ast_{OA} (EMD-CC)','X_{OA}(DL)','X_{OA}','X_{SIM} (SNR = 0.5)')
ylabel('Amplitude (\muV)')
ylim([-70 70])
xlim([1 Fs*time_length])
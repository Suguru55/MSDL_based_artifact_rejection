function [reconstructionsig_for_eog,reconstructionsig,...
          weights,icasig,Artifactual_IC_number] = automatic_wICA(input_data,Fs,ch,time)

%% ICA decomposition
% Infomax ICA
[weights,sphere,~,~,~,~,icasig] = runica(input_data);

%% Ocular artifact IC identification
% Calculate modified multiscale Fsle entropy(mMSE) 
% reference: H. B. Xie et al., 2008, "Measuring Time Series Regularity Using Nonlinear Similarity-Based Sample Entropy"

dim = 2;                                                                   % scale
freedom = ch;                                                              % degree of freedom

for channel = 1:ch
    data = icasig(channel,:);
    N = length(data(:,1:Fs*time));                                         % data reduction for fast calculation
    r = 0.2 * std(data(:,1:Fs*time),0,2);                                  % tolerance
    
    correl = zeros(1,2);
    dataMat = zeros(dim+1,N-dim);
    
    for i = 1:dim+1
        dataMat(i,:) = data(i:N-dim+i-1);
    end
    
    for m = dim:dim+1
        count = zeros(1,N-dim);
        tempMat = dataMat(1:m,:);
    
        for i = 1:N-m
            % calculate Chebyshev distance, excluding self-matching case
            dist = max(abs(tempMat(:,i+1:N-dim) - repmat(tempMat(:,i),1,N-dim-i)));
        
            %D = (dist < r);                                               % calculate Heaviside function of the distance
            D = 1./(1 + exp((dist-0.5)/r));                                % calculate Sigmoid function (mMSE)
            
            count(i) = sum(D)/(N-dim);
        end
        
        correl(m-dim+1) = sum(count)/(N-dim);
    end
    
    mMSE(channel,1) = log(correl(1)/correl(2));                           
end

% calculate kurtosis 
% references: A. Delorme et al., 2007, "Enhanced Detection of Artifacts in EEG Data Using Higher-Order Statistics and Independent Component Analysis"

Kurtosis = kurtosis(icasig(:,1:Fs*time)');
Kurtosis = Kurtosis';

% set boundaries
boundary_mMSE = mean(mMSE,1) - (std(mMSE,1)/sqrt(freedom)) * tinv(.95,freedom-1);
boundary_Kurtosis = mean(Kurtosis,1) + (std(Kurtosis,1)/sqrt(freedom)) * tinv(.95,freedom-1); 

% Identification
Artifactual_IC_number = [];

for i = 1:ch
    if (mMSE(i,:) <= boundary_mMSE) && (Kurtosis(i,:) >= boundary_Kurtosis)
        Artifactual_IC_number = [Artifactual_IC_number;i];
    end
end

%% Denoising by Discrete Wavelet Transform (DWT) 
% reference: N. P. Castellanos et al., 2006, "Recovering EEG Brain Signals: Artifact Suppression with Wavelet Enhanced Independent Component Analysis"

Denoised_icasig = icasig;                                                  % for eeg
Removed_icasig = zeros(ch,size(icasig,2));                                 % for artifact
Artifactual_IC = icasig(Artifactual_IC_number,:);

% Wavelet decomposition with Biorthogonal mother wavelet
Isdb4 = liftwave('db4');
els = {'p',[-0.125 0.125],0};
lsnew = addlift(Isdb4,els);

N = length(data);                  
                  
for i = 1:size(Artifactual_IC,1)
    [cA_1,cD_1] = lwt(Artifactual_IC(i,:),lsnew);         
    [cA_2,cD_2] = lwt(cA_1,lsnew);
    [cA_3,cD_3] = lwt(cA_2,lsnew);
    [cA_4,cD_4] = lwt(cA_3,lsnew);
    [cA_5,cD_5] = lwt(cA_4,lsnew);

    % \zeta(A) = sign(A)(abs(A)-T) : abs(A)>T   (soft)
    %              A               :            (hard)   
    %              0               : abs(A)<=T

    % Tj = median(|dj - median(dj)|)

    sigma1 = median(abs(cD_1))/0.6745;
    sigma2 = median(abs(cD_2))/0.6745;
    sigma3 = median(abs(cD_3))/0.6745;
    sigma4 = median(abs(cD_4))/0.6745;
    sigma5 = median(abs(cD_5))/0.6745;

    T1 = sigma1*sqrt(2*log(N));
    T2 = sigma2*sqrt(2*log(N));
    T3 = sigma3*sqrt(2*log(N));
    T4 = sigma4*sqrt(2*log(N));
    T5 = sigma5*sqrt(2*log(N));
    
    eeg_cD_1 = zeros(1,size(cD_1,2)); eeg_cD_2 = zeros(1,size(cD_2,2)); 
    eeg_cD_3 = zeros(1,size(cD_3,2)); eeg_cD_4 = zeros(1,size(cD_4,2));
    eeg_cD_5 = zeros(1,size(cD_5,2));
    artifact_cD_1 = zeros(1,size(cD_1,2)); artifact_cD_2 = zeros(1,size(cD_2,2));
    artifact_cD_3 = zeros(1,size(cD_3,2)); artifact_cD_4 = zeros(1,size(cD_4,2));
    artifact_cD_5 = zeros(1,size(cD_5,2));
    
    if max(abs(cD_1)) < T1
        %eeg_cD_1 = sign(cD_1).*(abs(cD_1) - T1);  % soft threshold
        eeg_cD_1 = cD_1;                         % hard threshold
    else
        artifact_cD_1 = cD_1;
    end
    
    if max(abs(cD_2)) < T2
        %eeg_cD_2 = sign(cD_2).*(abs(cD_2) - T2);
        eeg_cD_2 = cD_2;
    else
        artifact_cD_2 = cD_2;
    end
    
    if max(abs(cD_3)) < T3
        %eeg_cD_3 = sign(cD_3).*(abs(cD_3) - T3);
        eeg_cD_3 = cD_3;
    else
        artifact_cD_3 = cD_3;
    end
    
    if max(abs(cD_4)) < T4
        %eeg_cD_4 = sign(cD_4).*(abs(cD_4) - T4);
        eeg_cD_4 = cD_4;
    else
        artifact_cD_4 = cD_4;
    end
    
    if max(abs(cD_5)) < T5
        %eeg_cD_5 = sign(cD_5).*(abs(cD_5) - T5);
        eeg_cD_5 = cD_5;
    else
        artifact_cD_5 = cD_5;
    end
    
    eeg_cA_5 = zeros(1,size(cA_5,2));
    artifact_cA_5 = cA_5;
    
    % Inverse DWT
    [EEG_cA_4] = ilwt(eeg_cA_5,eeg_cD_5,lsnew);
    [EEG_cA_3] = ilwt(EEG_cA_4,eeg_cD_4,lsnew);
    [EEG_cA_2] = ilwt(EEG_cA_3,eeg_cD_3,lsnew);
    [EEG_cA_1] = ilwt(EEG_cA_2,eeg_cD_2,lsnew);
    [EEG_icasig] = ilwt(EEG_cA_1,eeg_cD_1,lsnew);
    
    [Artifact_cA_4] = ilwt(artifact_cA_5,artifact_cD_5,lsnew);
    [Artifact_cA_3] = ilwt(Artifact_cA_4,artifact_cD_4,lsnew);
    [Artifact_cA_2] = ilwt(Artifact_cA_3,artifact_cD_3,lsnew);
    [Artifact_cA_1] = ilwt(Artifact_cA_2,artifact_cD_2,lsnew);
    [Artifact_icasig] = ilwt(Artifact_cA_1,artifact_cD_1,lsnew);

    Denoised_icasig(Artifactual_IC_number(i,:),:) = EEG_icasig;
    Removed_icasig(Artifactual_IC_number(i,:),:) = Artifact_icasig;
end

%% Reconstruction
reconstructionsig = inv(weights * sphere) * Denoised_icasig; 
reconstructionsig_for_eog = inv(weights * sphere) * Removed_icasig;
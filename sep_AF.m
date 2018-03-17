%% Adaptive filter
%
% References: [1] P. He et al., 2004, "Removal of Ocular Artifacts from Electro-encephalogram
%                 by Adaptive Filtering"
%             [2] B. Yang et al., 2016 "Removal of EOG Artifacts from EEG Using a Cascade of 
%                 Sparse Autoencoder and Recursive Least Squares Adaptive Filter"
%             [3] https://jp.mathworks.com/matlabcentral/fileexchange/25769-adaptive-filter
%
% Nested function:
%
% Suguru Kanoga, last modification 22 Aug. 2017
% s.kanouga@aist.go.jp

function [results] = sep_AF(input_data,correct_eeg,testing_artifact,reference,Fs,time_length,code_dir)

M = 3;
lambda = 0.9999;
delta = 10^2;
ref_size = size(reference,1);

for sn_variation = 1:size(input_data,3)
    X = squeeze(input_data(:,:,sn_variation));
    [~,N] = size(X);    
    re_artifact_ep = zeros(size(X));
    re_eeg_ep = zeros(size(X));

    for ep = 1:N
        x = squeeze(reference(:,:,ep));  % reference inputs   
        d = X(:,ep);                     % primary input
        w = zeros(ref_size*M,1);
        R = diag((ones(ref_size*M,1))*delta);
        [i,j] = size(x);
             
        L_x = length(x);
        if (j>i)
            x = x.';
        end
        
        x = [zeros(M-1,ref_size); x];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % apply adaptive filter %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        for n = 2:L_x
            x_ref = x(n+M-1:-1:n,:); 
            x_ref = x_ref(:);
            
            R_n = R*x_ref;
            g = R_n/(lambda + x_ref'*R_n);
            Y = x_ref'*w;
            E = d(n,1) - Y;
           
            % updates
            w = w + g * E;
            R = (R - (g * x_ref' * R))/lambda;
            
            y(n,1) = x_ref'*w;
            e(n,1) = d(n,1) - y(n,1);
        end
        
        results.ProcessingTime(ep,sn_variation) = {toc};
        re_eeg_ep(:,ep) = e;
        re_artifact_ep(:,ep) = y;
        
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
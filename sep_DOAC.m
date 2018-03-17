%% Discrimitive ocular artifact correction
%
% reference: X. Li et al., 2017, "Discriminative Ocular Artifact Correction for Feature Learning in EEG Analysis"
%
% Nested functions:
%       - optimize_parameters.m (.../minFunc_2012/)
%       - minFunc               (.../minFunc_2012/minFunc/)
%       - mi                    (.../MIToolbox-3.0.0/matlab/)
%
% Suguru Kanoga, last modification 16 Aug. 2017
% s.kanouga@aist.go.jp

function [results] = sep_DOAC(input_data,correct_eeg,testing_artifact,Fs,time_length,code_dir,save_dir,class_tr)

% config
m = 10;
nbins = 12;
eta = 0.7;

for sn_variation = 1:size(input_data,3)
    temp = squeeze(input_data(:,:,sn_variation));
    [~,N] = size(temp);
    re_artifact_ep = zeros(size(temp));
    re_eeg_ep = temp;
    X_i_lib = zeros(3,Fs*time_length,N);
    
    tic
    for i = 1:N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % detect artfact segment in each epoch %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % moving average filter to obtain the smoothen signal x_s(t)
        x_s = movmean(temp(:,i),[m/2 m/2]);

        % calculate relative ampltude of the peaks
        % tau = 1;
        h = zeros(size(x_s));
        n_t = size(x_s,1);

        for t = 1:n_t
            if t == 1
                candidates = [0, abs(x_s(t,:)-x_s(t+1,:))];
            else
                if t == n_t
                    candidates = [abs(x_s(t,:)-x_s(t-1,:)), 0];
                else
                    candidates = [abs(x_s(t,:)-x_s(t-1,:)), 0, abs(x_s(t,:)-x_s(t+1,:))];
                end
            end
            
            h(t,:) = max(candidates);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define the peak amplitude rage parameter h_r %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obtain 12 peak amplitude bins
        % !! we modified the scale
        [p_c, edges] = histcounts(h,nbins);
        p_c = p_c./n_t;

        % h_b^1
        for bin = 1:size(p_c,2)
            if sum(p_c(:,1:bin)) > eta
                h_b1 = edges(:,1+bin);
                b1_bin = bin;
            
                break;
            end
        end
        
        % h_b^2 and h_u^1
        for bin = b1_bin+1:size(p_c,2)
            if sum(p_c(:,b1_bin:bin))/sum(p_c(:,b1_bin:end)) > eta
                h_b2 = edges(:,1+bin);
                h_u1 = h_b2;
                b2_bin = bin;
        
                break;
            end
        end
        
        % h_u^2
        for bin = b2_bin+1:size(p_c,2)
            if sum(p_c(:,b2_bin:bin))/sum(p_c(:,b2_bin:end)) > eta
                h_u2 = edges(:,1+bin);
        
                break;
            end
        end
        
        % find the set P_tj containing time indexes of those peaks with amplitude in
        % the range h_rj (= [h_bj h_uj]), j = [1,2].
        P_t1 = []; t_zb1 = []; t_za1 = [];
        P_t2 = []; t_zb2 = []; t_za2 = [];
    
        for t = 1:n_t
            if (m/2 < t && t < (n_t - m/2)) && (h_b1 < h(t,:) && h(t,:) < h_u1) 
                P_t1 = [P_t1; t];
            end
            
            if (m/2 < t && t < (n_t - m/2)) && (h_b2 < h(t,:) && h(t,:) < h_u2)
                P_t2 = [P_t2; t];
            end
        end
        
        % P_t.s have some cumbersome candidates, so let reduce the
        binary_x = x_s > 0;
        binary_x_delay = vertcat(binary_x(1),binary_x(1:end-1));
        diff_x = (binary_x_delay - binary_x) > 0;
        index = find(diff_x == 1); % zero-crossing points of x_s
    
        dammy = vertcat(0,P_t1(1:end-1));
        diff_t1 = (P_t1 - dammy) ~= 1;
        [row1,~] = find(diff_t1);
    
        for ind = 1:size(row1,1)
            index_b = index(find(index < P_t1(row1(ind),1)));
            index_a = index(find(index > P_t1(row1(ind),1)));
        
            t_zb1 = [t_zb1; index_b(knnsearch(index_b, P_t1(row1(ind),1)))];
            t_za1 = [t_za1; index_a(knnsearch(index_a, P_t1(row1(ind),1)))];
        end
        
        dammy = vertcat(0,P_t2(1:end-1));
        diff_t2 = (P_t2 - dammy) ~= 1;
        [row2,~] = find(diff_t2);
    
        for ind = 1:size(row2,1)
            index_b = index(find(index < P_t2(row2(ind),1)));
            index_a = index(find(index > P_t2(row2(ind),1)));
        
            t_zb2 = [t_zb2; index_b(knnsearch(index_b, P_t2(row2(ind),1)))];
            t_za2 = [t_za2; index_a(knnsearch(index_a, P_t2(row2(ind),1)))];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define X_a which contains all artifact signals x_ai %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % even if the time windows are overlapped, they will be linked together
        % after period extraction
        min_size = min(size(t_za1,1),size(t_zb1,1));
        window = [t_zb1(1:min_size), t_za1(1:min_size)];
    
        x_a1 = zeros(size(x_s));
        for ind = 1:size(window,1)
            x_a1(window(ind,1):window(ind,2)) = x_s(window(ind,1):window(ind,2));
        end
        
        min_size = min(size(t_za2,1),size(t_zb2,1));
        window = [t_zb2(1:min_size), t_za2(1:min_size)];
    
        x_a2 = zeros(size(x_s));
        for ind = 1:size(window,1)
            x_a2(window(ind,1):window(ind,2)) = x_s(window(ind,1):window(ind,2));
        end
        
        X_i = [temp(:,i), x_a1, x_a2]';
        X_i_lib(:,:,i) = X_i;
    end
    cd(save_dir);
    save('X_i_lib.mat', 'X_i_lib','class_tr');

    % change the order based on class labels
    X_i_1 = X_i_lib(:,:,class_tr==1);
    X_i_2 = X_i_lib(:,:,class_tr==2);
    X_i_lib_new = cat(3,X_i_1,X_i_2);
    [class_tr, ~] = sort(class_tr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimize the OAC coefficients %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % intialization
    theta_0= [1; 0; 0];
    lambda = [0 0.1 0.2 0.3 0.4 0.5];
    [a,b] = ndgrid(lambda,lambda);
    Lambda = [a(:),b(:)];
    Lambda = [Lambda, 1-b(:)]; % 36 combinations

    % Runs various limited-memory solvers on regularized_oscillatory_correction function for 25
    % function evaluations
    maxFunEvals = 25;
    options = [];
    options.display = 'none';
    options.maxFunEvals = maxFunEvals;
    
    for lam_num = 1:size(Lambda,1)
        global Lambda lam_num save_dir
        [theta_hat] = minFunc(@optimize_parameters,theta_0,options);
        if lam_num == 1
            [~, LASTID] = lastwarn;
            warning('off',LASTID);
        end
        
        % calculate I(f,c)
        f = zeros(N,3);
        for i = 1:N
            X_i = squeeze(X_i_lib_new(:,:,i));
            H_i = zeros(size(theta_0,1),n_t);
    
            for hil = 1:size(theta_0,1)
                H_i(hil,:) = imag(hilbert(X_i(hil,:)));
            end
            
            phi = sqrt((theta_hat' * X_i).^2 + (theta_hat' * H_i).^2);
            phi_bar = phi - trapz(phi)./n_t;
    
            % plus
            temp_i_plus = 0;
            for j = 1:sum(class_tr==1)
                X_plus = squeeze(X_i_lib_new(:,:,j));
                H_plus = zeros(size(theta_0,1),n_t);
                
                for hil = 1:size(theta_0,1)
                    H_plus(hil,:) = imag(hilbert(X_plus(hil,:)));
                end
                
                temp_i_plus =  temp_i_plus + sqrt((theta_hat' * X_plus).^2 + (theta_hat' * H_plus).^2);
            end
            
            phi_bar_plus = temp_i_plus./sum(class_tr==1);
            phi_bar_bar_plus = phi_bar_plus - trapz(phi_bar_plus)./n_t;
            
            %f(i,1) = (trapz(phi_bar.*phi_bar_bar_plus)/sqrt(trapz(phi_bar.^2).*trapz(phi_bar_bar_plus.^2)))./n_t;
            f(i,1) = trapz(phi_bar.*phi_bar_bar_plus)/sqrt(trapz(phi_bar.^2)*trapz(phi_bar_bar_plus.^2));
            
            % minus
            temp_i_minus = 0;
            for j = sum(class_tr==1)+1:sum(class_tr==1)+sum(class_tr==2)
                X_minus = squeeze(X_i_lib_new(:,:,j));
                H_minus = zeros(size(theta_0,1),n_t);
                
                for hil = 1:size(theta_0,1)
                    H_minus(hil,:) = imag(hilbert(X_minus(hil,:)));
                end
                
                temp_i_minus =  temp_i_minus + sqrt((theta_hat' * X_minus).^2 + (theta_hat' * H_minus).^2);
            end
        
            phi_bar_minus = temp_i_minus./sum(class_tr==2);
            phi_bar_bar_minus = phi_bar_minus - trapz(phi_bar_minus)./n_t;

            %f(i,2) = trapz(phi_bar.*phi_bar_minus)/sqrt(trapz(phi_bar.^2).*trapz(phi_bar_minus.^2))./n_t;
            f(i,2) = trapz(phi_bar.*phi_bar_bar_minus)/sqrt(trapz(phi_bar.^2).*trapz(phi_bar_bar_minus.^2));
            % power
            x_c = theta_hat' * X_i;
            f(i,3) = trapz(abs(x_c).^2) /n_t;
        end
        
        I = mi(f(:,1),class_tr) + mi(f(:,2),class_tr) + mi(f(:,3),class_tr);
    
        % normalize theta
        theta_hat = theta_hat./theta_hat(1);

        % initialize
        if lam_num == 1
            theta_0 = theta_hat;
            results.best_theta(:,sn_variation) = theta_hat;
            I_0 = I;
            results.best_lambda(sn_variation,:) = Lambda(lam_num,:);
            
%             disp('initial');
%             fprintf('I = %d\n',I);
%             disp(['theta = ',num2str(theta_hat')])
        else
            if sum(theta_0) == 1
                 results.best_theta(:,sn_variation) = theta_hat;
                 I_0 = I;
                 results.best_lambda(sn_variation,:) = Lambda(lam_num,:);
%                  disp('retake');
%                  fprintf('I = %d\n',I);
%                  disp(['theta = ',num2str(theta_hat')])
            end
        end
        
        % parameter update
        if (sum(theta_hat(:,1)) <= 1) && (theta_hat(2,1) >= -1) && (theta_hat(3,1) >= -1) ...
                && (theta_hat(2,1) <= 1) && (theta_hat(3,1) <= 1) && (I >= I_0) && (sum(theta_hat(2:3,1))~=0) % Eq. (25)   
            theta_0 = theta_hat;
            results.best_theta(:,sn_variation) = theta_hat;
            I_0 = I;
            results.best_lambda(sn_variation,:) = Lambda(lam_num,:);
            
%             disp('update');
%             fprintf('I = %d\n',I);
%             disp(['theta = ',num2str(theta_hat')])      
        end
    end
    results.LearningTime(sn_variation) = {toc};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reconstruct EEG and ocular signals %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:N
        tic
        artifact = results.best_theta(2,sn_variation)*squeeze(X_i_lib(2,:,i)) + results.best_theta(3,sn_variation)*squeeze(X_i_lib(3,:,i));
        eeg = temp(:,i)' + artifact;
    
        re_artifact_ep(:,i) = -artifact';
        re_eeg_ep(:,i) = eeg';
        results.ProcessingTime(i,sn_variation) = {toc};
        
        %%%%%%%%%%%%%%
        % evaluation %
        %%%%%%%%%%%%%%   
        results.mae_eeg(i,sn_variation) = {sum(abs(correct_eeg(:,i)-re_eeg_ep(:,i)))/(Fs*time_length)};
        results.mae_artifact(i,sn_variation) = {sum(abs(testing_artifact(:,i)-re_artifact_ep(:,i)))/(Fs*time_length)};
    end
    
    results.re_eeg(sn_variation) = {re_eeg_ep};
    results.re_artifact(sn_variation) = {re_artifact_ep};
end

cd(code_dir);
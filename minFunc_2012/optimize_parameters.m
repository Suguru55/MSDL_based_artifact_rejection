function [f,g] = optimize_parameters(theta)
global Lambda lam_num save_dir

cd(save_dir);
load('X_i_lib');
[~, n_t, N] = size(X_i_lib);

for cand = 1:size(theta,2)
    theta_sub = theta(:,cand);
   
    % optimaize parameters for each trial i
    for i = 1:N
        X_i = squeeze(X_i_lib(:,:,i));
        H_i = zeros(size(theta,1),n_t);
        
        for hil = 1:size(theta,1)
            H_i(hil,:) = imag(hilbert(X_i(hil,:)));
        end       
        
        phi = sqrt((theta_sub' * X_i).^2 + (theta_sub' * H_i).^2); % Eq. (13)
        phi_bar = phi - trapz(phi)./n_t; % Eq. (17)        
        
        % obtain an average oscilatory correlation between multiple trials
        j = 1:N; j(i) = [];
        oc = zeros(N-1,n_t);
    
        for l = 1:N-1
            X_j = X_i_lib(:,:,j(l));
            H_j = zeros(size(theta,1),n_t);
            
            for hil = 1:size(theta,1)
                H_j(hil,:) = imag(hilbert(X_j(hil,:)));
            end
            
            oc(l,:) = sqrt((theta_sub' *  X_j).^2 + (theta_sub' * H_j).^2);
        end
        
        psi = sum(oc)./(N-1);            % Eq. (14)
        psi_bar = psi - trapz(psi)./n_t; % Eq. (18)
    
        r_i_plus = 0; r_i_minus = 0;
        rho_hs_plus = 0; rho_hs_minus = 0;
    
        % plus: class 1, minus: class2
        if i < sum(class_tr==1)+1
            rho_hs_plus = rho_hs_plus + trapz(phi_bar.*psi_bar)/sqrt(trapz(phi_bar.^2)*trapz(psi_bar.^2)); % Eq. (16)
        
            for j = sum(class_tr==1)+1:sum(class_tr==1)+sum(class_tr==2)
                X_ic = squeeze(X_i_lib(:,:,j)); % for interclass
                H_ic = zeros(size(theta,1),n_t);
                
                for hil = 1:size(theta,1)
                    H_ic(hil,:) = imag(hilbert(X_ic(hil,:)));
                end

                phi_ic = sqrt((theta_sub' * X_ic).^2 + (theta_sub' * H_ic).^2);
                phi_bar_ic = phi_ic - trapz(phi_ic)./n_t;
                rho_hh = trapz(phi_bar.*phi_bar_ic)/sqrt(trapz(phi_bar.^2).*trapz(phi_bar_ic.^2));
                
                r_i_plus =  r_i_plus + rho_hh;
            end
        else
            rho_hs_minus = rho_hs_minus + trapz(phi_bar.*psi_bar)/sqrt(trapz(phi_bar.^2).*trapz(psi_bar.^2));
        
            for j = 1:sum(class_tr==1)
                X_ic = squeeze(X_i_lib(:,:,j));
                H_ic = zeros(size(theta,1),n_t);
                
                for hil = 1:size(theta,1)
                    H_ic(hil,:) = imag(hilbert(X_ic(hil,:)));
                end
                
                phi_ic = sqrt((theta_sub' * X_ic).^2 + (theta_sub' * H_ic).^2);
                phi_bar_ic = phi_ic - trapz(phi_ic)./n_t;
                rho_hh = trapz(phi_bar.*phi_bar_ic)/sqrt(trapz(phi_bar.^2).*trapz(phi_bar_ic.^2));
                
                r_i_minus = r_i_minus + rho_hh; 
            end
        end
        
        r_i(i,:) = r_i_plus + r_i_minus;
    end
    
    r_i = sum(r_i) / (sum(class_tr==1)*sum(class_tr==2)); % Eq. (19)

    % Eq. (20)
    r_w_plus = rho_hs_plus/sum(class_tr==1);
    r_w_minus = rho_hs_minus/sum(class_tr==2);

    % a regularized oscillatory correlation objective function
    f(cand) = (1-Lambda(lam_num,1))*(Lambda(lam_num,2) * r_w_plus + Lambda(lam_num,3) * r_w_minus) - Lambda(lam_num,1) * r_i; 
end

if nargout > 1 % gradiant required
    for cand = 1:size(theta,2)
        theta_size = size(theta,1);

        for var = 1:theta_size
            dth_r_i = zeros(N,1); % dth: dtheta
            for i = 1:N
                theta_sub = theta(:,cand);
                dth_X_i = squeeze(X_i_lib(var,:,i));
                dth_H_i = imag(hilbert(dth_X_i));
           
                dth_phi = (theta_sub(var) .* (dth_X_i.^2) + theta_sub(var) .* (dth_H_i.^2)) ./ ...
                    sqrt((theta_sub(var) .* dth_X_i).^2 + (theta_sub(var) .* dth_H_i).^2);
                dth_phi(isnan(dth_phi)) = 0;
                dth_phi_bar = dth_phi - trapz(dth_phi)./n_t; 
                
                % obtain an average oscilatory correlation between multiple trials
                j = 1:N;
                j(i) = [];
                dth_oc = zeros(N-1,n_t);
    
                for l = 1:N-1
                    dth_X_j = X_i_lib(var,:,j(l));
                    dth_H_j = imag(hilbert(dth_X_j));
                    
                    dth_oc(l,:) = (theta_sub(var) .* (dth_X_j.^2) + theta_sub(var) .* (dth_H_j.^2)) ./ ...
                                     sqrt((theta_sub(var) .* dth_X_j).^2 + (theta_sub(var) * dth_H_j).^2);
                end
                dth_oc(isnan(dth_oc)) = 0;
                
                dth_psi = sum(dth_oc)./(N-1);
                dth_psi_bar = dth_psi - trapz(dth_psi)./n_t;
    
                dth_r_i_plus = 0; dth_r_i_minus = 0;
                dth_rho_hs_plus = 0; dth_rho_hs_minus = 0;
    
                if i < sum(class_tr==1)+1
                    dth_rho_hs_plus = dth_rho_hs_plus + trapz(dth_phi_bar.*dth_psi_bar)/sqrt(trapz(dth_phi_bar.^2).*trapz(dth_psi_bar.^2));
                    dth_rho_hs_plus(isnan(dth_rho_hs_plus)) = 0;
        
                    for j = sum(class_tr==1)+1:sum(class_tr==1)+sum(class_tr==2)
                        dth_X_ic = squeeze(X_i_lib(var,:,j)); % for interclass
                        dth_H_ic = imag(hilbert(dth_X_ic));
                        
                        dth_phi_ic = (theta_sub(var) .* (dth_X_ic.^2) + theta_sub(var) .* (dth_H_ic.^2)) ./ ...
                            sqrt((theta_sub(var) .* dth_X_ic).^2 + (theta_sub(var) .* dth_H_ic).^2);
                        dth_phi_bar_ic = dth_phi_ic - trapz(dth_phi_ic)./n_t;
                        dth_rho_hh = trapz(dth_phi_bar.*dth_phi_bar_ic)/sqrt(trapz(dth_phi_bar.^2).*trapz(dth_phi_bar_ic.^2));
                        dth_rho_hh(isnan(dth_rho_hh)) = 0;
                        
                        dth_r_i_plus =  dth_r_i_plus + dth_rho_hh;
                    end
                else
                    dth_rho_hs_minus = dth_rho_hs_minus + trapz(dth_phi_bar.*dth_psi_bar)/sqrt(trapz(dth_phi_bar.^2).*trapz(dth_psi_bar.^2));
                    dth_rho_hs_minus(isnan(dth_rho_hs_minus)) = 0;
                    
                    for j = 1:sum(class_tr==1)
                        dth_X_ic = squeeze(X_i_lib(var,:,j)); % for interclass
                        dth_H_ic = imag(hilbert(dth_X_ic));
                        
                        dth_phi_ic = (theta_sub(var) .* (dth_X_ic.^2) + theta_sub(var) .* (dth_H_ic.^2)) ./ ...
                            sqrt((theta_sub(var,:) .* dth_X_ic).^2 + (theta_sub(var,:) .* dth_H_ic).^2);
                        dth_phi_bar_ic = dth_phi_ic - trapz(dth_phi_ic)./n_t;
                        dth_rho_hh = trapz(dth_phi_bar.*dth_phi_bar_ic)/sqrt(trapz(dth_phi_bar.^2).*trapz(dth_phi_bar_ic.^2));
                        dth_rho_hh(isnan(dth_rho_hh)) = 0;
                        
                        dth_r_i_minus = dth_r_i_minus + dth_rho_hh; 
                    end
                end
                
                dth_r_i(i,:) = dth_r_i_plus + dth_r_i_minus;
            end
            
            dth_r_i = sum(dth_r_i) / (sum(class_tr==1)*sum(class_tr==2));

            dth_r_w_plus = dth_rho_hs_plus/sum(class_tr==1);
            dth_r_w_minus = dth_rho_hs_minus/sum(class_tr==2);

            % a regularized oscillatory correlation objective function
            g(var,cand) = (1-Lambda(lam_num,1))*(Lambda(lam_num,2) * dth_r_w_plus + Lambda(lam_num,3) * dth_r_w_minus) - Lambda(lam_num,1) * dth_r_i; 
            
            if g(var,cand) == 0
               g(var,cand) = 1; 
            end
        end
    end
end
%% SingleChannel EMD with Cross Correlation
%
% reference: R. Patel et al., 2016, "Suppression of Eye-blink Associated Artifact Using Single
%            Channel EEG Data by Combining Cross-correlatiuon with
%            Empirical Mode Decomposition"
%
% Nested function:
% - emd
%
% !!Identification step was arraged.
% Suguru Kanoga, last modification 16 Aug. 2017
% s.kanouga@aist.go.jp

function [results] = sep_EMDandCC(input_data,correct_eeg,testing_artifact,Fs,time_length,code_dir)

window_length = 1.5; % 1.5 s data segment
MaxIter = 5000;
threshold = 0.5;

for sn_variation = 1:size(input_data,3)
    X = squeeze(input_data(:,:,sn_variation));
    [~,N] = size(X);
    freedom = ceil(time_length/window_length);
    re_artifact_ep = zeros(size(X));
    re_eeg_ep = X;

    for ep = 1:N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % detect artifact segment in each epoch %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        artifact_seg = [];
        sub_pos = [];
        seg_lib = zeros(Fs*window_length,ceil(time_length/window_length));
        tic
    
        for win = 1:ceil(time_length/window_length)
            if win == ceil(time_length/window_length)
                seg_X = X(1+(win-1)*Fs*window_length:end,ep);
            else
                seg_X = X(1+(win-1)*Fs*window_length:win*Fs*window_length,ep);
            end
            
            temp_v(win,:) = var(seg_X);
            temp_s(win,:) = skewness(seg_X);
        
            seg_lib(1:size(seg_X,1),win) = seg_X;
        end
        
        % !! all epochs have already been identified as artifactual by visual
        %    inspection from BCI competition
        %    ¨ each epoch must contain at least one artufactual segment
        %      but the artifact segment (1.5 s) may present once or twice, sometimes three times in
        %      the epoch (6 s)
        %      ¨ set boundaries based on t-value (the DOF is ceil(time_length/window_length) - 1)
        boundary_v = mean(temp_v,1) - std(temp_v,1)/sqrt(freedom) * tinv(.95,freedom-1);
        boundary_s = mean(temp_s,1) - std(temp_s,1)/sqrt(freedom) * tinv(.95,freedom-1);
    
        for i = 1:ceil(time_length/window_length)
            if (temp_v(i,:) >= boundary_v) && (temp_s(i,:) >= boundary_s)
                artifact_seg = [artifact_seg; i];
            
                % get position infomation that has maximum value in artifactual segment
                [~,temp] = max(seg_lib(:,i));
                sub_pos = [sub_pos; temp];
            end
        end
        
        pos = sub_pos + (Fs*window_length)*(artifact_seg-1);
    
        % move the position to 200 ms (50-th sampling point) and save templates
        % ¨ if a position value is larger than (Fs*time_length - Fs + (Fs/10)*3) or smaller than Fs/5,
        %   the segmenet should be zero-padded

        for j = 1:size(artifact_seg,1)
            ex_seg = [];
            if pos(j,1) >= Fs*time_length - (Fs+(Fs/10)*3)
                sub_ex_seg = X(pos(j,1)-(Fs/5-1):end,ep);
                ex_seg = [sub_ex_seg; zeros(Fs+(Fs/2) - size(sub_ex_seg,1),1)];
            else
                if pos(j,1) <= Fs/5
                    sub_ex_seg = X(1:pos(j,1)+(Fs+(Fs/10)*3),ep);
                    ex_seg = [zeros(Fs+(Fs/2) - size(sub_ex_seg,1),1); sub_ex_seg];
                else
                    ex_seg = X(pos(j,1)-(Fs/5-1):pos(j,1)+(Fs+(Fs/10)*3),ep);
                end
            end
            template = ex_seg(1+Fs/10:(Fs/10)*3,:);
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % apply EMD for the segments %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            desvio_estandar = std(ex_seg);
            ex_seg = ex_seg / desvio_estandar;

            [modos, ~, ~] = emd(ex_seg,'MAXITERATIONS',MaxIter);
            modos = modos * desvio_estandar;                               % IMFs
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculate CC between template and IMFs %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            modos_art = zeros(size(modos));
        
            for k = 1:size(modos,1)
                cc = xcorr(modos(k,:),template)/(norm(modos(k,:))*norm(template));
                if max(cc) >= threshold      
                    modos_art(k,:) = modos(k,:);
                    modos(k,:) = zeros;
                end
            end
            
            % reconstruction
            re_seg_artifact = sum(modos_art);
            re_seg_eeg = sum(modos);
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % realign reconstructed segments into original segments %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if pos(j,1) >= Fs*time_length - (Fs+(Fs/10)*3)
                re_eeg_ep(pos(j,1)-(Fs/5-1):end,ep) = re_seg_eeg(1,1:size(sub_ex_seg,1));
                re_artifact_ep(pos(j,1)-(Fs/5-1):end,ep) = re_seg_artifact(1,1:size(sub_ex_seg,1));
            else
                if pos(j,1) <= Fs/5
                    re_eeg_ep(1:pos(j,1)+(Fs+(Fs/10)*3),ep) = re_seg_eeg(1,1+Fs+(Fs/2) - size(sub_ex_seg,1):end);
                    re_artifact_ep(1:pos(j,1)+(Fs+(Fs/10)*3),ep) = re_seg_artifact(1,1+Fs+(Fs/2) - size(sub_ex_seg,1):end);               
                else
                    re_eeg_ep(pos(j,1)-(Fs/5-1):pos(j,1)+(Fs+(Fs/10)*3),ep) = re_seg_eeg;
                    re_artifact_ep(pos(j,1)-(Fs/5-1):pos(j,1)+(Fs+(Fs/10)*3),ep) = re_seg_artifact;
                end
            end
        end
        results.ProcessingTime(ep,sn_variation) = {toc};
        
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
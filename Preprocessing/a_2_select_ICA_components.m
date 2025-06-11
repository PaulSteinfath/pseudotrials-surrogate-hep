function a_2_select_ICA_components(data_path, save_path, error_path, qa_path, ica_window, thresholds, markers)
    % A_2_SELECT_ICA_COMPONENTS Remove ICA components based on various artifact criteria
    %
    % The function:
    %   1. Detects R-peaks for ECG component identification
    %   2. Identifies and removes ICA components related to:
    %      - Heart artifacts (ECG)
    %      - Muscle artifacts
    %      - Eye movements
    %      - Line noise
    %      - Channel noise
    %      - Other artifacts
    %   3. Generates quality assessment plots
    %   4. Saves cleaned datasets and error logs
    %
    % Inputs:
    %   data_path   - Path to input .set files.
    %   save_path   - Path to save cleaned datasets.
    %   error_path  - Path to save error logs.
    %   qa_path     - Path to save quality assessment plots.
    %   ica_window  - Window size for ICA analysis.
    %   thresholds  - Structure containing thresholds for artifact rejection:
    %                 * thresholds('sd_ecg') - ECG component rejection threshold.
    %                 * thresholds('ecg')    - ECG component correlation threshold.
    %                 * thresholds('muscle') - Muscle artifact rejection threshold.
    %                 * thresholds('eye')    - Eye movement artifact rejection threshold.
    %                 * thresholds('ln')     - Line noise rejection threshold.
    %                 * thresholds('chann')  - Channel noise rejection threshold.
    %                 * thresholds('other')  - Other artifact rejection threshold.
    %   markers     - Structure containing marker names:
    %                 * markers('r') - R-peak marker name.
    %
    % Outputs:
    %   No direct outputs. Saves to disk:
    %   - Cleaned EEG datasets (.set files) in save_path
    %   - Quality assessment plots in qa_path/[subject_id]/
    %     * PSD pre/post ICA
    %     * Removed components topoplots
    %     * All components topoplots
    %     * Potential brain components topoplots
    %     * ECG-correlated components plots
    %   - Error logs in error_path (if errors occur)
    
    %% get file names
    files = dir(fullfile(data_path, '*.set'));
    
    %% get already processed data
    savefiles = dir(fullfile(save_path, '*.set'));
    savefilesnames = {savefiles.name};
    
    %% Loop over subjects
    for i=1:length(files)
        try
            subjid = files(i).name(1:end-4);
    
            % if subject exists, skip
            if any(strcmp(savefilesnames, [subjid, '.set']))
                continue
            end
            
            % create QA folders
            if not(isfolder([qa_path, subjid]))
                mkdir([qa_path, subjid])
            end
    
            % load data
            EEG = pop_loadset('filename',files(i).name, 'filepath', data_path);
    
            %% Detect R-Peaks
            % Detect R-peaks for clarity in identifying ECG components
            [EEG, ~]  = insert_R_peaks(EEG, EEG.ECG.data(1,:), 1, markers('r'));
    
            % keep original, non-epoched EEG
            EEGorig = EEG;
    
            % Select ICA components based on ECG template
            [EEG, cV, rejV] = ecg_ica_corr(EEG, [], ica_window, markers('r'), thresholds('sd_ecg'));
            
            % plot & save correlated ECG components
            plot_ecg_ica_comps(EEG, rejV,cV, [qa_path, subjid, '/'], subjid);
            close all
            
            % run IClabel for classification of components on continuous data
            EEGorig = pop_iclabel(EEGorig, 'default');
    
            %% Heart Components 
            % Use correlation + ICLabel to find components
            % correlated with the ECG > ecg_tresh and not classified as brain by IClabel
            rej_heart_pos =  find(rejV + (cV>thresholds('ecg')) == 2); 
            rej_heart_pos_corr = rej_heart_pos(~(EEGorig.etc.ic_classification.ICLabel.classifications(rej_heart_pos,1)>0.99));
            EEGorig.etc.ic_remove.heart = rej_heart_pos_corr; 
    
            %% Muscle, eye, line noise, channel noise and other Components
            % Add above threshold component idx to original set
            EEGorig.etc.ic_remove.muscle = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,2) >= thresholds('muscle'));
            EEGorig.etc.ic_remove.eye = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,3) >= thresholds('eye'));
            EEGorig.etc.ic_remove.LineNoise = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,5) >= thresholds('ln'));
            EEGorig.etc.ic_remove.ChannelNoise = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,6) >= thresholds('chann'));      
            EEGorig.etc.ic_remove.other = find(EEGorig.etc.ic_classification.ICLabel.classifications(:,7) >= thresholds('other'));
            
            % save correlation and index of potential heart components
            EEGorig.etc.ECG.ic_ecg_corr = cV;
            EEGorig.etc.ECG.ic_ecg_corr_above_thres = rejV;
            
            % reject components
            reject = unique([EEGorig.etc.ic_remove.heart', EEGorig.etc.ic_remove.muscle', EEGorig.etc.ic_remove.eye', EEGorig.etc.ic_remove.LineNoise', EEGorig.etc.ic_remove.ChannelNoise', EEGorig.etc.ic_remove.other']);
            EEGclean = pop_subcomp(EEGorig, reject);
            EEGclean = eeg_checkset(EEGclean);
           
            % save cleaned dataset
            EEGclean = pop_saveset(EEGclean,'filename',subjid,'filepath',save_path);
    
            %% Some quality assessment (QA) Plots
            
            % Plot spectogram
            figure; plot_spec(EEGclean.data(1:31,:)', EEG.srate, 'f_max',70);
            saveas(gcf, [qa_path, subjid, '/', subjid, '_PSD_post_ICA.png'],'png');
            
            % pre-clean PSD
            figure; plot_spec(EEG.data(1:31,:)', EEG.srate, 'f_max',70);
            saveas(gcf, [qa_path, subjid, '/', subjid, '_PSD_pre_ICA.png'],'png');
    
            % plot removed components
            pop_topoplot(EEGorig, 0, reject ,subjid,[] ,0,'electrodes','on','iclabel', 'on', 'iclabel', 'on');
            saveas(gcf, [qa_path, subjid, '/', subjid, '_removed_components.png'],'png');
    
            % plot all components
            pop_topoplot(EEGorig, 0, [1:size(EEGorig.icaact)] ,subjid,[] ,0,'electrodes','on', 'iclabel', 'on');
            saveas(gcf, [qa_path, subjid, '/', subjid, '_all_components.png'],'png');
    
            % plot artefact components that have high probability to be brain activity
            brain_comps = reject(EEGorig.etc.ic_classification.ICLabel.classifications(reject,1)>0.995);
            if any(brain_comps)
                pop_topoplot(EEGorig, 0,brain_comps,subjid,[] ,0,'electrodes','on', 'iclabel', 'on');
                saveas(gcf, [qa_path, subjid, '/', subjid, '_potential_brain_components.png'],'png');
            end
            
            % if Loop completed and previously a error log was created, delete the error log
            if isfile([error_path, subjid, '_error_log_ICA_comp_select.txt'])
                delete([error_path, subjid, '_error_log_ICA_comp_select.txt'])
            end
                
        catch ME
    
            fileID = fopen([error_path, subjid, '_error_log_ICA_comp_select.txt'],'w');
            fprintf(fileID,'%6s\n',subjid);  
            fprintf(fileID,'%6s\n',getReport(ME));
            fclose(fileID);          
            
        end
            
    end
end

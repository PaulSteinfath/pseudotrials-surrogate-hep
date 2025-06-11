function a_3b_apply_ICA_components_to_filtered_data_add_task_marker(task_raw_path, rest_raw_path, save_path, rest_pre_ica_path, rest_post_ica_path, error_path, qa_path, fs, elecfile, highpass_cu, lowpass_cu, line_noise_f)
    % a_3b_apply_ICA_components_to_filtered_data_add_task_marker - Apply
    % identified ICA components to data filtered in a different frequency
    % range and add markers from task to resting state data for analysis
    %
    % This function applies identified Independent Component Analysis (ICA) components
    % to EEG data that has been filtered in a different frequency range. It also adds
    % task markers from task data to resting state data for further analysis.
    % It processes each subject's data by loading raw and ICA data, applying
    % filters, merging events, interpolating channels, removing noisy ICA components,
    % and saving the results. It also generates quality assurance plots for each subject.
    %
    % Inputs:
    %   - task_raw_path: Path to the raw task data files.
    %   - rest_raw_path: Path(s) to the raw resting state data files.
    %   - save_path: Path where the processed data will be saved.
    %   - rest_pre_ica_path: Path to the pre-ICA processed data files.
    %   - rest_post_ica_path: Path to the post-ICA processed data files.
    %   - error_path: Path where error logs will be saved.
    %   - qa_path: Path where quality assurance outputs will be saved.
    %   - fs: Sampling frequency for resampling the data.
    %   - elecfile: File containing electrode information for channel editing.
    %   - highpass_cu: High-pass filter cutoff frequency.
    %   - lowpass_cu: Low-pass filter cutoff frequency.
    %   - line_noise_f: Frequency of line noise to be notched out.
    %
    % Outputs:
    %   - Processed EEG data files saved to the specified save_path.
    %   - Quality assurance plots saved to the qa_path.
    %   - Error logs saved to the error_path if any issues occur during processing.
    %

    
    %% Get files
    
    % Raw oddball data - for task markers
    filesOddball = dir(fullfile(task_raw_path, '*.vhdr'));
    filesOddball = {filesOddball.name};
    
    % Raw resting state data - to filter and preprocess again
    rawFiles = dir(fullfile(rest_raw_path, '*.vhdr'));
    rawFilesList = {rawFiles.name};
    
    prepFiles = dir(fullfile(save_path, '*.set'));
    prepFileList = {prepFiles.name};
    
    % Files containing ICA matrices
    ICAFiles = dir(fullfile(rest_pre_ica_path, '*.set'));
    ICAFilesList = {ICAFiles.name};
    
    % Files containing info about which components to reject
    postICAFiles = dir(fullfile(rest_post_ica_path, '*.set'));
    postICAFilesList = {postICAFiles.name};
    
    %% Loop over subjects
%     parfor s=1:length(postICAFiles)
    for s=1:length(postICAFiles)
        try 
            % get subject iD
            subjid = postICAFiles(s).name(1:end-4);
            newID = [subjid '.set'];     
            rawID = [subjid '.vhdr'];
            oddballID = postICAFiles(s).name(1:10);
    
            % Check if subject was processed already, if so skip
            if any(strcmp(prepFileList, newID))
                continue
            end
    
            %% match post ICA and Raw data
            post_ica_ID =  find(strcmp(postICAFilesList, newID));
            ica_ID =  find(strcmp(ICAFilesList, newID));
            raw_ID = find(strcmp(rawFilesList, rawID));
            raw_ID_oddball = find(strcmp(filesOddball, [oddballID, '_Novelty.vhdr']));
    
            % If subjids dont match - give error
            if ~(strcmp(postICAFilesList{post_ica_ID}(1:end-4), rawFiles(raw_ID).name(1:end-5)) && strcmp(rawFiles(raw_ID).name(1:end-5), ICAFilesList{ica_ID}(1:end-4)))   
                fileID = fopen([error_path, subjid, 'no_match_IDs_error_log.txt'],'w');
                fprintf(fileID,'%6s\n',subjid, postICAFilesList{post_ica_ID}(1:end-4), ICAFilesList{ica_ID}(1:end-4));
                fprintf(fileID,'%6s\n',getReport(ME));
                fclose(fileID);
                continue
            end
            
            %% load all the data
            %load ICA data
            EEG_ica = pop_loadset(ICAFilesList{ica_ID}, rest_pre_ica_path);
            
            %load post ICA data
            EEG_post_ica = pop_loadset(postICAFilesList{post_ica_ID}, rest_post_ica_path);
            
            %load Raw
            [EEG, com] = pop_loadbv(rest_raw_path, rawFiles(raw_ID).name);
    
            %load Oddball
            [EEG_oddball, com] = pop_loadbv(task_raw_path, filesOddball{raw_ID_oddball});
    
            %% resample to fs
            EEG = pop_resample(EEG,fs); 
            EEG_oddball = pop_resample(EEG_oddball,fs); 
    
            %% merge the events
            EEG.event = [EEG_oddball.event, EEG.event];
            eeg_checkset(EEG);
    
            %% do the preprocessing again on raw data
            EEG.setname = rawFiles(raw_ID).name;
    
            % Add montage 
            EEG=pop_chanedit(EEG, 'lookup',elecfile);
             
            % Band-pass filter and notch (using standardized cutoffs)
            [b,a] = butter(2, highpass_cu/(EEG.srate/2), 'high');
            EEG.data = filtfilt(b,a, double(EEG.data)')';
            [c,d] = butter(2, lowpass_cu/(EEG.srate/2), 'low');
            EEG.data = filtfilt(c,d, double(EEG.data)')';
            
            % Notch filter
            nf_lowpass_cu = line_noise_f - 1;
            nf_highpass_cu = line_noise_f + 1;
            EEG = pop_eegfiltnew(EEG, 'locutoff', nf_lowpass_cu, 'hicutoff', nf_highpass_cu, 'revfilt', 1, 'plotfreqz', 0);
            
            % resample to fs
            EEG = pop_resample(EEG, fs);
            
            % Select channels
            EEG.ECG = pop_select(EEG, 'channel', {'EKG'});
            EEG.VEOG = pop_select(EEG, 'channel', {'VEOG'});
            EEG.HEOG = pop_select(EEG, 'channel', {'HEOG'});
            EEG = pop_select(EEG, 'channel', [1:31]);
    
            %keep backup for channel interpolation
            originalEEG=EEG;
            
            %% remove flat channels and interpolate
            EEG = pop_select(EEG, 'nochannel',EEG_ica.reject.removed_channels);
    
            allchan = {originalEEG.chanlocs.labels};
            EEG.reject.removed_channels = allchan(~ismember({originalEEG.chanlocs.labels}, {EEG.chanlocs.labels}));
            
            %interpolate
            EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
            eeg_checkset(EEG);   
            
            %% Add ICA components and remove noisy ones
            EEG.icawinv = EEG_ica.icawinv;
            EEG.icasphere = EEG_ica.icasphere;
            EEG.icaweights = EEG_ica.icaweights;
            EEG.icachansind = EEG_ica.icachansind;         
            eeg_checkset(EEG);        
            
            reject = unique([EEG_post_ica.etc.ic_remove.heart', EEG_post_ica.etc.ic_remove.muscle', EEG_post_ica.etc.ic_remove.eye', EEG_post_ica.etc.ic_remove.LineNoise', EEG_post_ica.etc.ic_remove.ChannelNoise', EEG_post_ica.etc.ic_remove.other']);
            EEG = pop_subcomp(EEG, reject);
            EEG.etc = EEG_post_ica.etc;
            EEG = eeg_checkset(EEG);
                                 
            %% save results
            pop_saveset(EEG,'filename',subjid,'filepath',save_path);
          
            %% QA steps
            
            % Plot spectogram
            figure; plot_spec(EEG.data', fs, 'f_max',70);
            saveas(gcf, [qa_path, subjid, '/', subjid, '_PSD_post_03_45Hz_ICA.png'],'png');
            
            %if Loop completed and previously a error log was created, delete it
            if isfile([error_path, subjid, '_error_log.txt'])
                delete([error_path, subjid, '_error_log.txt'])
            end
            
        catch ME   
            fileID = fopen([error_path, subjid, '_error_log.txt'],'w');
            fprintf(fileID,'%6s\n',subjid);  
            fprintf(fileID,'%6s\n',getReport(ME));
            fclose(fileID);          
        end
    end
end

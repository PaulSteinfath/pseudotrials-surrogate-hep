function a_1_oddball_preprocessing(data_path, save_path, error_path, fs, elecfile, ica_highpass_cu, ica_lowpass_cu,line_noise_f, flatline_crit)
    % A_1_ODDBALL_PREPROCESSING - Preprocess EEG data for oddball experiment
    %
    % This function performs the initial preprocessing steps for EEG data:
    %   - Band-pass filtering (high-pass and low-pass)
    %   - Notch filtering for line noise
    %   - Downsampling to a specified sampling rate
    %   - Channel selection and removal of bad channels
    %   - Artifact removal and interpolation
    %   - Independent Component Analysis (ICA)
    %
    % Inputs:
    %   data_path        - Path to directory containing raw EEG files (string)
    %   save_path         - Path to directory for saving preprocessed EEG files (string)
    %   error_path        - Path to directory containing error logs (string)
    %   fs               - Target sampling rate for downsampling (double)
    %   elecfile         - Path to electrode montage file (string)
    %   ica_highpass_cu  - High-pass filter cutoff for ICA preprocessing (double)
    %   ica_lowpass_cu   - Low-pass filter cutoff for ICA preprocessing (double)
    %   line_noise_f     - Line noise frequency for notch filtering (e.g., 50 Hz or 60 Hz) (double)
    %   flatline_crit    - Criterion for detecting flatline channels (double)
    %
    % Outputs:
    %   No direct outputs. Preprocessed EEG data is saved to `savepath` in `.set` format.

    %% get files
    files = dir(fullfile(data_path, '*.vhdr'));

    %% get (if any) already preprocessed filenames 
    prepFiles = dir(fullfile(save_path, '*.set'));
    prepFileList = {prepFiles.name};

    %% preprocess data and save as set
    parfor i = 1:length(files)
    
        try 
                 
            % subject ID and name for saving
            subjid = files(i).name(1:end-5);
            newID = [subjid '.set'];
            
    
            % if subject is already in prep-folder, skip it
            if any(strcmp(prepFileList, newID))
                 continue
            end
    
            % load data
            [EEG, com] = pop_loadbv(data_path, files(i).name);
            
            % keep subjid stored in EEG struct
            EEG.setname = files(i).name(1:end-5);
    
            % add montage EEG electrodes
            EEG=pop_chanedit(EEG, 'lookup',elecfile);
         
    
            % Band-pass filter combined with the notch
            [b,a]=butter(2,ica_highpass_cu/(EEG.srate/2),'high'); % high pass filter
            EEG.data=filtfilt(b,a,double(EEG.data)')'; 
    
            [c,d]=butter(2,ica_lowpass_cu/(EEG.srate/2), 'low'); % low pass filter
            EEG.data=filtfilt(c,d,double(EEG.data)')';
            
            % Notch filter for line-noise filtering
            nf_lowpass_cu = line_noise_f - 1;
            nf_highpass_cu = line_noise_f + 1;
            EEG = pop_eegfiltnew(EEG, 'locutoff',nf_lowpass_cu,'hicutoff',nf_highpass_cu,'revfilt',1,'plotfreqz',0);
    
            % resample to fs
            EEG = pop_resample(EEG, fs); 
    
            % select ECG and EOG channels and save it in its own struct
            % inside EEG struct
            EEG.ECG =   pop_select(EEG, 'channel',{'EKG'});
            EEG.VEOG =  pop_select(EEG, 'rmchannel',{'VEOG'}); 
            EEG.HEOG =  pop_select(EEG, 'rmchannel',{'HEOG'}); 
            EEG =       pop_select(EEG, 'channel',[1:31]);
    
    
            % keep backup for channel interpolation
            originalEEG=EEG;
            
            % remove flat channels
            EEG = clean_artifacts(EEG,'ChannelCriterion','off', 'FlatlineCriterion', flatline_crit, 'BurstCriterion', 'off', 'WindowCriterion','off');
            
            % store removed channels in EEG struct
            allchan = {originalEEG.chanlocs.labels};
            EEG.reject.removed_channels = allchan(~ismember({originalEEG.chanlocs.labels}, {EEG.chanlocs.labels}));
    
            % interpolate
            EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
            eeg_checkset(EEG);        
                  
            % get rank of the data
            dataRank = sum(eig(cov(double(EEG.data'))) > 1E-7);
    
            % run extended infomax ICA
            EEG = pop_runica(EEG, 'icatype', 'runica', 'pca', dataRank , 'options', {'extended' 1});
                
            % save results
            pop_saveset(EEG,'filename',subjid,'filepath',save_path);
    
        catch ME           
            % in case something goes wrong, save an error log
            fileID = fopen([error_path, subjid, '_error_log.txt'],'w');
            fprintf(fileID,'%6s\n',subjid);  
            fprintf(fileID,'%6s\n',getReport(ME));
            fclose(fileID);
                    
        end
    end
end

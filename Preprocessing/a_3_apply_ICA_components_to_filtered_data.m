function a_3_apply_ICA_components_to_filtered_data(raw_path, save_path, pre_ica_path, post_ica_path, error_path, qa_path, fs, elecfile, highpass_cu, lowpass_cu, line_noise_f, chan_crit, ln_crit)
    % A_3_APPLY_ICA_COMPONENTS_TO_FILTERED_DATA - Applies ICA components to filtered EEG data.
    %
    % Description:
    %    This function processes EEG data by applying Independent Component Analysis (ICA)
    %    components to originally filtered data. It performs the following steps:
    %    1) Loads raw, pre-ICA, and post-ICA EEG data files.
    %    2) Matches and verifies subject IDs across different data files.
    %    3) Resamples, filters, and interpolates EEG data.
    %    4) Applies ICA weights and removes identified artifact components.
    %    5) Cleans the data using the clean_rawdata plugin.
    %    6) Saves the processed EEG data and generates QA plots.
    %    7) Handles errors and logs them if necessary.
    %
    % Inputs:
    %    raw_path - Path to the raw EEG data files (.vhdr format).
    %    save_path - Path to save the filtered EEG data files.
    %    pre_ica_path - Path to the pre-ICA processed EEG data files.
    %    post_ica_path - Path to the post-ICA processed EEG data files.
    %    error_path - Path to save error logs.
    %    qa_path - Path for saving quality assurance (QA) plots.
    %    fs - Sampling frequency for resampling the EEG data.
    %    elecfile - File containing electrode montage information.
    %    highpass_cu - High-pass filter cutoff frequency.
    %    lowpass_cu - Low-pass filter cutoff frequency.
    %    line_noise_f - Line noise frequency for notch filtering (e.g., 50 Hz or 60 Hz) (double)
    %    chan_crit - Channel criterion for clean_rawdata plugin.
    %    ln_crit - Line noise criterion for clean_rawdata plugin.
    %
    % Outputs:
    %    No direct outputs. The function saves filtered EEG data files and QA plots to specified directories.


    
    %% get all the file names
    rawFiles = dir(fullfile(raw_path, '*.vhdr')); 
    prepFiles = dir(fullfile(save_path, '*.set')); 
    icaFiles = dir(fullfile(pre_ica_path, '*.set'));
    postICAFiles = dir(fullfile(post_ica_path, '*.set'));
    
    rawFilesList = {rawFiles.name};
    prepFileList = {prepFiles.name};
    icaFilesList = {icaFiles.name};
    postICAFilesList = {postICAFiles.name};

    %% loop over subjects
    for i=1:length(postICAFilesList)
        try 
                    
            % get subject iD
            subjid = postICAFiles(i).name(1:end-4);
            newID = [subjid '.set'];     
            rawID = [subjid '.vhdr'];          
            
            % Check if subject was processed already - if so skip it            
            if any(strcmp(prepFileList, newID))
                continue
            end
    
            % match post ICA and Raw data
            post_ica_ID =  find(strcmp(postICAFilesList, newID));
            ica_ID =  find(strcmp(icaFilesList, newID));
            raw_ID = find(strcmp(rawFilesList, rawID));
            
            % If subjids dont match - give error
            if ~(strcmp(postICAFilesList{post_ica_ID}(1:end-4), rawFiles(raw_ID).name(1:end-5)) && strcmp(rawFiles(raw_ID).name(1:end-5), icaFilesList{ica_ID}(1:end-4)))   
                fileID = fopen([save_error, subjid, 'no_match_IDs_error_log.txt'],'w');
                fprintf(fileID,'%6s\n',subjid, postICAFilesList{post_ica_ID}(1:end-4), icaFilesList{ica_ID}(1:end-4));
                fprintf(fileID,'%6s\n',getReport(ME));
                fclose(fileID);
                continue
            end
    
            % load ICA data
            EEG_ica = pop_loadset(icaFilesList{ica_ID}, pre_ica_path);
            
            % load post ICA data
            EEG_post_ica = pop_loadset(postICAFilesList{post_ica_ID}, post_ica_path);
            
            % load Raw
            [EEG, com] = pop_loadbv(raw_path, rawFiles(raw_ID).name);
    
            % resample to sampling_rate
            EEG = pop_resample(EEG, fs); 
    
            EEG.setname = rawFiles(i).name(1:end-5);
    
            % Add montage 
            EEG=pop_chanedit(EEG, 'lookup',elecfile);
             
            % select channels
            EEG.ECG =   pop_select(EEG, 'channel',{'EKG'});
            EEG.VEOG =  pop_select(EEG, 'rmchannel',{'VEOG'}); 
            EEG.HEOG =  pop_select(EEG, 'rmchannel',{'HEOG'}); 
            EEG =       pop_select(EEG, 'channel',[1:31]);
    
            % remove & interpolate flat channels
            EEG = pop_select(EEG, 'nochannel',EEG_ica.reject.removed_channels);
            EEG = pop_interp(EEG, EEG_ica.chanlocs, 'spherical');
            
            % Keep the info
            EEG.reject.removed_channels = EEG_ica.reject.removed_channels;
                    
            % low pass than high pass filter the data
            [b,a]=butter(2,highpass_cu/(EEG.srate/2),'high'); % high pass filter
            EEG.data=filtfilt(b,a,double(EEG.data)')'; 
            [c,d]=butter(2,lowpass_cu/(EEG.srate/2), 'low'); % low pass filter
            EEG.data=filtfilt(c,d,double(EEG.data)')';
            
            % Notch filter for line-noise filtering
            nf_lowpass_cu = line_noise_f - 1;
            nf_highpass_cu = line_noise_f + 1;
            EEG = pop_eegfiltnew(EEG, 'locutoff',nf_lowpass_cu,'hicutoff',nf_highpass_cu,'revfilt',1,'plotfreqz',0);
    
            % keep backup for channel interpolation
            originalEEG=EEG;
                    
            % Add ICA results to EEG struct
            EEG.icawinv = EEG_ica.icawinv;
            EEG.icasphere = EEG_ica.icasphere;
            EEG.icaweights = EEG_ica.icaweights;
            EEG.icachansind = EEG_ica.icachansind;         
            eeg_checkset(EEG);        
                   
            % remove ICA components
            reject = unique([EEG_post_ica.etc.ic_remove.heart', EEG_post_ica.etc.ic_remove.muscle', EEG_post_ica.etc.ic_remove.eye', EEG_post_ica.etc.ic_remove.LineNoise', EEG_post_ica.etc.ic_remove.ChannelNoise', EEG_post_ica.etc.ic_remove.other']);
            EEG = pop_subcomp(EEG, reject);
            EEG.etc = EEG_post_ica.etc;
            EEG = eeg_checkset(EEG);
    
            % remove noisy & low correlated channels using clean_rawdata plugin
            EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion',chan_crit,'LineNoiseCriterion',ln_crit,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
    
            allchan = {originalEEG.chanlocs.labels};
            EEG.reject.removed_channels = allchan(~ismember({originalEEG.chanlocs.labels}, {EEG.chanlocs.labels}));
    
            % interpolate
            EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
            eeg_checkset(EEG);      
            
            % save results
            pop_saveset(EEG,'filename',subjid,'filepath',save_path);
            
            %% QA
            % Plot spectogram
            figure; plot_spec(EEG.data', fs, 'f_max',70);
            saveas(gcf, [qa_path, subjid, '/', subjid, '_PSD_post_', highpass_cu, '_', lowpass_cu, 'Hz_ICA.png'],'png');
                    
            close all
            
            % If Loop completed and previously a error log was created, delete the error log
            if isfile([error_path, subjid, '_error_log.txt'])
                delete([error_path, subjid, '_error_log.txt'])
            end
                    
        catch ME
            % in case something goes wrong, save an error log
            fileID = fopen([error_path, subjid, '_error_log.txt'],'w');
            fprintf(fileID,'%6s\n',subjid);  
            fprintf(fileID,'%6s\n',getReport(ME));
            fclose(fileID);
                        
        end
    end
end

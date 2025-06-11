function a_4_sort_random_pseudo(data_path, save_path, error_path, qa_path, binlister_file, epoch_length,lda_file, ids_file, markers, min_isi, isi, hep_prestim_window, hep_epoch, noise_epoch_tresh) 
    % a_4_sort_random_pseudo: Process EEG data for pseudo-random sorting
    %
    % This function processes EEG data files by applying various preprocessing
    % steps, including de-meaning, event checking, R-peak detection, reaction
    % time calculation, and more. It also performs median split ERP image
    % analysis using LDA and saves the processed data in ERPLab format.
    %
    % Inputs:
    % - data_path: Path to the directory containing EEG data files (.set).
    % - save_path: Path to save the processed ERP files.
    % - error_path: Path to save error logs if any processing fails.
    % - qa_path: Path for quality assurance output directories.
    % - binlister_file: File containing binlister information for ERPLab.
    % - epoch_length: Length of the epoch for processing.
    % - lda_file: File containing LDA matrices.
    % - ids_file: File containing subject IDs.
    % - markers: Struct containing event markers.
    % - min_isi: Minimum inter-stimulus interval.
    % - isi: Inter-stimulus interval.
    % - hep_prestim_window: Pre-stimulus window for HEP analysis.
    % - hep_epoch: HEP epoch length.
    % - noise_epoch_tresh: Noise threshold for epoch rejection.
    %
    % Outputs:
    % - Processed ERP files saved in the specified save_path directory.
    % - Error logs saved in the error_path directory if processing fails.
    % - Quality assurance output directories created in the qa_path for each subject.
    
    prepFileList = dir(fullfile(save_path, '*.erp'));
    
    %% Prep settings
    
    % subjects to process
    files = dir(fullfile(data_path, '*.set'));
    
    % Load LDA matrices
    erp_t_filter = struct2array(load(lda_file));
    ids = loadtxt(ids_file);
    
    %% Loop over subjects
    % parfor i=1:length(files)
    for i=1:length(files)
        try
    
            % get subjid
            subjid = files(i).name(1:end-4);
            
            % check if file already exists, if so skip it
            if any(strcmp({prepFileList.name}, [subjid '.erp']))
                continue
            end
    
            % create qa output dirs
            create_dirs({qa_path}, subjid)
    
            %% load data
            EEG = pop_loadset('filename',files(i).name,'filepath', data_path);
    
            %% de-mean the signal
            EEG = pop_rmbase( EEG, [],[]);
    
            %% Check if task events are present
            EEG = check_events(EEG, {markers('standard_stim'), markers('target_stim'), markers('novelty_stim')}, markers('response'), min_isi, error_path, subjid);
    
            %% Detect R-Peaks and add them to EEG
            [EEG, ~] = insert_R_peaks(EEG, EEG.ECG.data(1,:), 1, markers('r'));
    
            %% Add reaction time
            EEG = calcRT(EEG, markers('target_stim'), markers('response'), isi);
            
            EEG = matchHEPtoTask_diffPseudo(EEG, markers('target_stim'), hep_prestim_window, markers('r'), markers('pseudo'), markers('random')); % 'S 10', 'S 20', 'S 30'}
            
            %% split by RT
            EEG = split_by_RT_pseudoHEP(EEG, markers);
            EEG = split_by_RT_target(EEG, markers); % Note: If I understand well, this is not used (already done in split_by_RT_pseudoHEP with different markers) and could be removed
            % Check events consistency
            EEG = eeg_checkset(EEG, 'eventconsistency');
    
            %% split by activity
            % Match LDA id's with EEG id's
            LDA_idx = find(contains(ids, subjid(1:10)));
            filter = erp_t_filter(LDA_idx, :)* 10e-7; % rescale filter weights to microvolt
         
            % Define combined markers for specific conditions
            % marker for heartbeat before target (8820)
            r_before_target_marker = [markers('r') markers('target')];
            % marker for target with heartbeat in the pre-stim time (8822)
            target_with_r_prestim_marker = [markers('r') markers('before_target')];
            % marker for pseudo heartbeat before target (48820)
            pseudo_r_before_target_marker = [markers('pseudo') markers('r') markers('target')];
            % marker for target with pseudo heartbeat (48822)
            target_with_pseudo_r_prestim_marker = [markers('pseudo') markers('r') markers('before_target')];
            % marker for random heartbeat before target (88820)
            random_r_before_target_marker = [markers('random') markers('r') markers('target')];
            % marker for target with random heartbeat (88822)
            target_with_random_r_prestim_marker =[markers('random') markers('r') markers('before_target')];
    
            event_pairs = {
                {target_with_r_prestim_marker, r_before_target_marker}, ...
                {target_with_pseudo_r_prestim_marker, pseudo_r_before_target_marker}, ...
                {target_with_random_r_prestim_marker, random_r_before_target_marker},
                };
           
            % using LDA
            EEG = median_split_erpimage_fake_HEP_random_pseudo(EEG, hep_epoch, epoch_length, filter, event_pairs, binlister_file, markers);
    
            % check event consistency
            EEG = eeg_checkset(EEG,'eventconsistency');
    
            %% prepare data for epoching in ERPLab format, epoch the data
            EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' } ); % GUI: 23-Oct-2020 11:02:15
            EEG  = pop_binlister( EEG , 'BDF', binlister_file, 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG', 'UpdateEEG', 'on' ); % GUI: 23-Oct-2020 15:00:27
    
            % No baseline  correction -  separate for Task HEP
            EEG = pop_epochbin(EEG , [epoch_length(1)  epoch_length(2)], 'none'); % GUI: 23-Oct-2020 15:02:15
    
            %% clean epochs, remove non-zero events, remove duplicates & empty epochs
            EEG = remove_non_zero_events(EEG);
    
            %% Remove Noise
            EEG  = pop_artmwppth( EEG , 'Channel',  1:31, 'Flag',  1, 'Threshold',  noise_epoch_tresh, 'Twindow', [epoch_length(1)  epoch_length(2)], 'Windowsize',  200, 'Windowstep',  100 ); % GUI: 25-Jan-2022 14:13:20
            EEG = pop_rejepoch( EEG, find(EEG.reject.rejmanual),0);
         
            %% Create average ERP set
            ERP = pop_averager( EEG , 'Criterion', 'good', 'DQ_flag', 1, 'DSindex', 1, 'ExcludeBoundary', 'on', 'SEM', 'on' );% GUI: 23-Oct-2020 15:20:27
    
            ERP = pop_savemyerp(ERP, 'erpname', subjid, 'filename', [subjid, '.erp'], 'filepath', save_path);
            
            %% if Loop completed and previously a error log was created, delete the error logs
            if isfile([error_path, subjid, '_HEP_error_log.txt'])
                delete([error_path, subjid, '_HEP_error_log.txt']);
            end
    
        catch ME
            fileID = fopen([error_path, subjid, '_HEP_error_log.txt'],'w');
            fprintf(fileID,'%6s\n',subjid);
            fprintf(fileID,'%6s\n',getReport(ME));
            fclose(fileID);
        end
    
    end
end

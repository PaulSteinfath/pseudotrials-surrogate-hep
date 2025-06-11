function [EEG, HEP_qrs] = insert_R_peaks(EEG, ECG, add_ecg, r_marker)
    %% Insert R-Peak Markers into EEG Data
    % This function adds R-peak markers to EEG data by optionally appending ECG data to the EEG dataset.
    % It detects R-peaks using the `heplab_slowdetect` function and saves these events into the EEG structure.
    %
    % Inputs:
    %   EEG      - EEG data structure
    %   ECG      - ECG data to be added to the EEG
    %   add_ecg  - Boolean flag indicating whether to add ECG data to EEG (1 to add, 0 to skip)
    %
    % Outputs:
    %   EEG      - Updated EEG data structure with R-peak events
    %   HEP_qrs  - Detected R-peak events
    
    % add back ECG to EEG
    if add_ecg == 1
        % Append ECG data as a new channel in the EEG data
        EEG.data(end+1,:,:) = ECG;
        % Update the number of channels in the EEG structure
        EEG.nbchan = size(EEG.data,1);
        % Label the new channel as 'ECG'
        EEG.chanlocs(end+1).labels = 'ECG';
        % check consistency of fields in EEG
        EEG = eeg_checkset(EEG);
    end
    
    % Detect R-peaks in ECG data
    HEP_qrs = heplab_slowdetect(ECG,EEG.srate);
    
    % Save detected R-peak events to EEG set
    EEG = fun_heplab_save_events(EEG, HEP_qrs, r_marker);
    % check consistency of events in EEG
    EEG = eeg_checkset(EEG,'eventconsistency');

end
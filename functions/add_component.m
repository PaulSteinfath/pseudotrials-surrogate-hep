function EEG = add_component(EEG, event_idx, timing, amplitude, srate, npoints)
    % ADD_COMPONENT Adds a waveform component to EEG data at a specified event.
    %
    %   EEG = ADD_COMPONENT(EEG, EVENT_IDX, TIMING, AMPLITUDE, SRATE, NPOINTS)
    %   adds a waveform component to the EEG data structure at the event
    %   specified by EVENT_IDX. The waveform is defined by its TIMING (a two-element
    %   vector specifying onset and duration in seconds), AMPLITUDE, sampling rate
    %   SRATE, and the total number of data points NPOINTS.
    %
    %   Inputs:
    %   - EEG: EEG data structure containing event information.
    %   - EVENT_IDX: Index of the event in the EEG.event structure.
    %   - TIMING: Two-element vector [onset, duration] in seconds.
    %   - AMPLITUDE: Amplitude of the waveform to be added.
    %   - SRATE: Sampling rate of the EEG data.
    %   - NPOINTS: Total number of data points in the EEG data.
    %
    %   Outputs:
    %   - EEG: Updated EEG data structure with the added waveform component.
    %
    %   The function generates a waveform using a Gaussian shape modulated by a
    %   Hann window, and adds it to the EEG data at the specified event location.
    %   The onset and duration are rounded to ensure valid indices, and the function
    %   ensures that the waveform is added within the bounds of the data.

    % Get base latency (floating point in EEGLAB)
    base_latency = EEG.event(event_idx).latency;
    
    % Round to ensure valid indices
    onset = round(base_latency + (timing(1) * srate));
    duration = round(timing(2) * srate);
    
    % Enforce positive onset and valid end index
    onset = max(1, onset);
    end_point = min(npoints, onset + duration - 1);
    
    % Ensure onset < end_point
    if onset >= end_point
        return; % Skip adding if there's no valid range
    end
    
    % Generate waveform (example Gaussian * Hann window)
    waveform_length = end_point - onset + 1;
    gauss_center = waveform_length / 2;
    gauss_sd = waveform_length / 4;
    
    % Time vector from 1 to waveform_length
    t = 1:waveform_length;
    
    % Gaussian shape
    gauss_wave = exp(-((t - gauss_center).^2) / (2 * gauss_sd^2));
    
    % Hann window
    hann_window = hann(waveform_length)';
    
    % Combine and scale
    waveform = amplitude * gauss_wave .* hann_window;
    
    % Add to EEG data
    EEG.data(:, onset:end_point) = EEG.data(:, onset:end_point) + waveform;
end
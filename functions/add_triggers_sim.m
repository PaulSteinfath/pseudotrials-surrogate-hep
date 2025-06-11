function EEG = add_triggers_sim(EEG, stim_interval, hep_prestim_window, stim_trigger, r_trigger)
    % ADD_TRIGGERS_SIM Adds simulated triggers to EEG data for testing purposes.
    %
    % Description:
    %   This function adds two types of triggers to EEG data:
    %   1. Stimulus triggers at regular intervals
    %   2. Heartbeat triggers before every second stimulus with random delays
    %
    % Inputs:
    %   EEG               - EEGLAB EEG structure containing the data
    %   stim_interval     - Interval between stimuli in seconds
    %   hep_prestim_window - Two-element vector specifying the range [min, max] in seconds for random delays before stimulus
    %   stim_trigger      - String representing the stimulus trigger type (e.g., 'S 20')
    %   r_trigger         - String representing the heartbeat trigger type (e.g., '88')
    %
    % Outputs:
    %   EEG               - Modified EEG structure with added triggers
    %
    % Trigger Types:
    %   stim_trigger - Stimulus trigger (e.g., 'S 20')
    %   r_trigger    - Heartbeat trigger (e.g., '88')
    %
    % Notes:
    %   - Heartbeat triggers are added only before even-numbered stimuli
    %   - Heartbeat triggers are placed randomly between the specified range before stimulus
    %   - The function requires the EEGLAB toolbox for eeg_checkset
    %
    % Example:
    %   EEG = add_triggers_sim(EEG, 2, [0.6, 1.1], 'S 20', '88'); % Add triggers with 2-second intervals
    
    %% Initialize parameters
    srate = EEG.srate;
    stimulus_interval = srate * stim_interval;  % Convert interval to data points
    event_counter = 1;
    events = struct('type', {}, 'latency', {}, 'urevent', {}); 
    
    %% Add triggers throughout the recording
    trial_counter = 0;
    
    for t = stimulus_interval:stimulus_interval:size(EEG.data,2)
        trial_counter = trial_counter + 1;
        
        % Add heartbeat trigger before even-numbered stimuli
        if mod(trial_counter, 2) == 0
            % Calculate random delay inside hep_prestim_window
            pre_stim_delay = rand() * hep_prestim_window(1);
            % Ensure trigger is at least 600ms before stimulus and within data bounds
            pre_stim_onset = max(1, t - round(pre_stim_delay * srate) - round(hep_prestim_window(2) * srate));
            
            % Add heartbeat event
            events(event_counter) = struct(...
                'type', r_trigger, ...
                'latency', pre_stim_onset, ...
                'urevent', event_counter);
            event_counter = event_counter + 1;
        end
        
        % Add stimulus event
        events(event_counter) = struct(...
            'type', stim_trigger, ...
            'latency', t, ...
            'urevent', event_counter);
        event_counter = event_counter + 1;
    end
    
    %% Update EEG structure
    EEG.event = events;
    EEG = eeg_checkset(EEG); % Validate the modified dataset
end

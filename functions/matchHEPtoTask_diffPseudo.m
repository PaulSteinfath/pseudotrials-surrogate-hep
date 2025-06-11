function EEG = matchHEPtoTask_diffPseudo(EEG, task_events, hep_time_window, r_marker, pseudo_marker, random_marker)
    %   Matches heart evoked potentials (HEP) to task events and simulates
    %   pseudo-events. 
    %   The function loops through each specified task event, finds the closest preceding heartbeat within the
    %   specified time window, and updates the event markers accordingly. If no heartbeat is found, a simulated
    %   event is created. Random triggers are also added for each event as pseudotrials.  
    %
    %   Inputs:
    %       EEG - EEG data structure containing event information.
    %       task_events - Cell array of event types to match with heartbeats.
    %       hep_time_window - Two-element vector specifying the time window (in seconds) before each event to search for heartbeats.
    %       r_marker - String specifying the type of heartbeat event marker.
    %       pseudo_marker - String specifying simulated event.
    %       random_marker - String specifying random event.
    %
    %   Outputs:
    %       EEG - Updated EEG data structure with matched and simulated events.
    %

    % Define markers
    % pseudo heartbeat marker
    pseudo_r_marker = [pseudo_marker r_marker];
    % random trigger marker
    random_trig_marker = [random_marker r_marker];
    
    % Loop through each specified event type
    for evnt = 1:length(task_events)
        event = task_events{evnt};
    
        % Find indices of the specific event
        idx_event = find(contains({EEG.event.type}, event));
        lat_event = cell2mat({EEG.event(idx_event).latency});
    
        % Find indices of heartbeat events (r_marker)
        idx_r_marker = find(contains({EEG.event.type}, r_marker));
        lat_r_marker = cell2mat({EEG.event(idx_r_marker).latency});
    
        % Initialize counters for matched and random triggers
        num_matched = 0;
        num_random = 0;
    
        % Loop through each occurrence of the event
        for ii = 1:length(idx_event)
            hep_start = lat_event(ii) - EEG.srate * hep_time_window(1);
            hep_end = lat_event(ii) - EEG.srate * hep_time_window(2);
    
            % Find heartbeats within the specified time window
            fit_heps = find(lat_r_marker > hep_start & lat_r_marker < hep_end);
    
            % Select the last heartbeat before the event
            if length(fit_heps) > 1
                fit_heps = fit_heps(end);
            end
    
            % Ensure no other heartbeats occurred between the selected heartbeat and the event
            if ~isempty(fit_heps) && any(lat_r_marker(fit_heps) < lat_r_marker & lat_r_marker < lat_event(ii))
                fit_heps = [];
            end
    
            % Handle events based on their type and whether a suitable heartbeat was found
            if isempty(fit_heps)
                % No suitable heartbeat found; simulate one
                rand_lat = hep_start + (hep_end - hep_start) * rand;
                EEG.event(end+1) = EEG.event(idx_r_marker(1));
                EEG.event(end).latency = rand_lat;
                EEG.event(end).type = [pseudo_r_marker, event(3:end)];
    
                EEG.event(end+1) = EEG.event(idx_event(ii));
                EEG.event(idx_event(ii)).type = [pseudo_r_marker, event(3:end-1) '2'];
    
                num_random = num_random + 1;
            else
                % Suitable heartbeat found; update types accordingly
                EEG.event(end+1) = EEG.event(idx_r_marker(fit_heps));
                EEG.event(idx_r_marker(fit_heps)).type = [r_marker, event(3:end)];
    
                EEG.event(end+1) = EEG.event(idx_event(ii));
                EEG.event(idx_event(ii)).type = [r_marker, event(3:end-1) '2'];
    
                num_matched = num_matched + 1;
            end
    
        % Always add a random trigger as an alternative pseudotrial
        rand_lat = hep_start + (hep_end - hep_start) * rand;
        EEG.event(end+1) = EEG.event(idx_r_marker(1));
        EEG.event(end).latency = rand_lat;
        EEG.event(end).type = [random_trig_marker, event(3:end)];
    
        EEG.event(end+1) = EEG.event(idx_event(ii));
        EEG.event(idx_event(ii)).type = [random_trig_marker, event(3:end-1) '2'];
    
        end
    
    end
end

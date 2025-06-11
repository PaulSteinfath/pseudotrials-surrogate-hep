function [EEG, reaction_times] = calcRT(EEG, target_stim, response_stim, isi)
% calcRT - Calculates reaction time based on target and response stimuli.
%
% INPUTS:
%   EEG           - EEGLAB data structure.
%   target_stim   - Name of the target stimulus (e.g., 'S 20').
%   response_stim - Name of the response stimulus (e.g., 'S 1').
%   isi           - Inter-stimulus interval (time allowed for participant to respond in seconds).
%
% OUTPUTS:
%   EEG           - EEGLAB data structure with reaction times added.
%   reaction_times - Array of calculated reaction times in milliseconds.

% Find indices and latencies for response and target stimuli
idx_button_press = find(strcmp({EEG.event.type}, response_stim));
lat_button_press = [EEG.event(idx_button_press).latency];

idx_target = find(strcmp({EEG.event.type}, target_stim));
lat_target = [EEG.event(idx_target).latency];

reaction_times = [];

% Calculate reaction times
for ii = 1:length(lat_target)
    
    % Determine the time window for valid responses
    response_time_limit = lat_target(ii) + EEG.srate * isi;
    valid_responses = lat_button_press(lat_button_press > lat_target(ii) & lat_button_press < response_time_limit);

    % If no valid response is found, continue to the next target
    if isempty(valid_responses)
        continue;
    end
    
    % Use the first valid response (closest to the target stimulus)
    first_response_latency = valid_responses(1);
    reaction_time = first_response_latency - lat_target(ii);
    
    % Store the reaction time in milliseconds
    reaction_time_ms = (reaction_time / EEG.srate) * 1000;
    reaction_times = [reaction_times, reaction_time_ms];
    
    % Add reaction time to the EEG.event structure for both response and target
    EEG.event(idx_button_press(lat_button_press == first_response_latency)).RT = reaction_time_ms;
    EEG.event(idx_target(ii)).RT = reaction_time_ms;
end

% Check event consistency
EEG = eeg_checkset(EEG, 'eventconsistency');

end

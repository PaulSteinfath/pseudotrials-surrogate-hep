function [EEG, err_msg] = check_events(EEG, task_events, responses, min_isi, save_error, subjid)
% CHECK_EVENTS - Validates the presence of specified events in EEGLAB data
% structure and checks inter-stimulus intervals (ISI).
%
% INPUTS:
% EEG        - EEGLAB data structure.
% task_events - Cell array of task-related events to be checked.
% responses  - Cell array of response-related events to be checked.
% min_isi    - Minimum inter-stimulus interval allowed (in seconds).
% save_error - Path where error logs should be saved.
% subjid     - Subject identifier for logging purposes.
%
% OUTPUTS:
% EEG        - Modified EEGLAB data structure with events removed if ISI is violated.
% err_msg    - Error message(s) generated during processing.

% Initialize empty error message
err_msg = [];

% Flatten nested cell arrays if necessary
if iscell(task_events{1})
    task_events = [task_events{:}];
end
if iscell(responses{1})
    responses = [responses{:}];
end

% Combine all events to check (task events and responses)
all_events = [task_events, responses];
ispresent = cellfun(@(s) any(strcmp({EEG.event.type}, s)), all_events);

% Check if necessary events are present
missing_events = all_events(~ispresent);
if ~isempty(missing_events)
    err_msg = sprintf('No %s events were found.', strjoin(missing_events, ', '));
    log_error(save_error, subjid, '_no_task_events_error_log.txt', err_msg);
end

% Find indices of task events
task_event_idx = find(ismember({EEG.event.type}, task_events));

% Check inter-stimulus interval (ISI) and remove events that violate the ISI
if ~isempty(task_event_idx)
    isi = diff([EEG.event(task_event_idx).latency] / EEG.srate);
    bad_marker_pos = task_event_idx(isi < min_isi);

    if ~isempty(bad_marker_pos)
        err_msg = sprintf('A total of %d event marker(s) were removed due to ISI < %.2f sec.', length(bad_marker_pos), min_isi);
        EEG.event(bad_marker_pos) = [];  % Remove events with ISI violations
        log_error(save_error, subjid, '_events_removed_error_log.txt', err_msg);
    end
end

end

function log_error(save_path, subjid, log_name, message)
% LOG_ERROR - Helper function to log error messages.
%
% INPUTS:
% save_path - Path to save the log file.
% subjid    - Subject identifier for logging purposes.
% log_name  - Name of the log file.
% message   - Error message to log.

fileID = fopen(fullfile(save_path, [subjid, log_name]), 'w');
fprintf(fileID, 'Subject: %s\nError: %s\n', subjid, message);
fclose(fileID);
end
function EEG = remove_non_zero_events(EEG)
    % REMOVE_NON_ZERO_EVENTS - Removes non-zero latency events from EEG data.
    %
    % Input:
    %   EEG - EEG structure containing epoch and event information.
    %
    % Output:
    %   EEG - EEG structure with only zero latency events retained.
    %
    % Description:
    %   This function processes the EEG structure to retain only events that
    %   occur at latency 0 within each epoch. It removes duplicate events
    %   and ensures event consistency. Additionally, it removes empty epochs
    %   from the EEG data.
    %
    % Notes:
    %   - The function checks for the presence of specific event-related fields
    %     and updates them accordingly.
    %   - It uses 'eeg_checkset' to ensure event consistency after modifications.
    %   - Empty epochs are identified and removed using 'pop_select'.
    %
    % See also: eeg_checkset, pop_select

    %keep only events at latency 0, remove all else since there are many duplicates
    for kk = 1:size(EEG.epoch,2)
           
        if isfield(EEG, 'epoch') && numel(EEG.epoch(kk).event) > 1
        
            keep_evnt_idx = cell2mat(EEG.epoch(kk).eventlatency) == 0;
        
            event_fields = {'event', 'eventbepoch', 'eventbini', 'eventbinlabel', 'eventbvmknum', 'eventbvtime', ...
                            'eventchannel', 'eventcodelabel', 'eventduration', 'eventenable', 'eventflag', ...
                            'eventitem', 'eventlatency', 'eventtype', 'eventvisible'};
        
            for i = 1:length(event_fields)
                field_name = event_fields{i};
                if isfield(EEG.epoch(kk), field_name)
                    EEG.epoch(kk).(field_name) = EEG.epoch(kk).(field_name)(keep_evnt_idx);
                end
            end        
        end    
    end
  
    %events to keep
    EEG.event = EEG.event([EEG.epoch.event]);
    EEG = eeg_checkset(EEG,  'eventconsistency' );
    
    % remove duplicate events (since events are multiplied during epoching, remove the arbitrary ones and keep one event per epoch)
    if isfield(EEG.event, 'item')
        [~, ia, ~] = unique([EEG.event.item]);
        EEG.event = EEG.event(ia);
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end
    
    % remove empty epochs
    emptyEpochs = find(cellfun(@isempty,{EEG.epoch.event}));
    EEG = pop_select(EEG, 'notrial', emptyEpochs);


end
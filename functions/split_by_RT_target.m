function EEG = split_by_RT_target(EEG, markers)
    % Split according to RT
    % Extract specific markers from the input 'markers' structure
    % target stim marker
    target_stim_marker = markers('target_stim');
    % target marker
    target_marker = markers('target');
    % marker for trials with low amplitude target
    low_target_marker = markers('low_target');
    % marker for trials with high amplitude target
    high_target_marker = markers('high_target');

    % Define combined markers for specific conditions
    target_low_rt_marker = [target_marker low_target_marker];
    target_high_rt_marker = [target_marker high_target_marker];


    % Find indices of events matching specific markers
    idx_erp = find(ismember({EEG.event.type},target_stim_marker));
    RTs = [EEG.event(idx_erp).RT];
    
    med_sort_val = median(RTs);
    low_RT = find(RTs<med_sort_val);
    high_RT = find(RTs>med_sort_val);
    
    EEG_event_temp = EEG.event; %backup of EEG.event
    [EEG_event_temp(idx_erp(low_RT)).type] = deal(target_low_rt_marker);
    [EEG_event_temp(idx_erp(high_RT)).type] = deal(target_high_rt_marker);
    
    keep_indices = [idx_erp(low_RT), idx_erp(high_RT)];
    EEG.event = cat(2,EEG.event,EEG_event_temp(keep_indices)); %append low/high erps to original event structure, so that the original events are preserved


end
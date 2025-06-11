function EEG = split_by_RT_pseudoHEP(EEG, markers)
    % SPLIT_BY_RT_PSEUDOHEP Splits EEG events by reaction time and categorizes them.
    %
    % This function processes EEG event data by splitting events into low and high
    % reaction time (RT) categories based on the median RT. It updates the event
    % types for Heartbeat Evoked Potentials (HEP) and Event-Related Potentials (ERP)
    % trials, including pseudo and random pseudo conditions.
    % It identifies events based on predefined markers, splits them into
    % low and high RT categories, and updates the event types accordingly. It
    % processes three types of trials: HEP, pseudo HEP, and random pseudo HEP.
    %
    % INPUTS:
    %   EEG     - A structure containing EEG data, including an 'event' field with
    %             event types and reaction times.
    %   markers - A function or structure that returns specific marker codes when
    %             given a marker name as input.
    %
    % OUTPUT:
    %   EEG     - The input EEG structure with updated event types reflecting the
    %             categorization into low and high RT events for different trial
    %             conditions.
    %

    %% Define markers
    % Extract specific markers from the input 'markers' structure
    % R marker
    r_marker = markers('r');
    % pseudo marker
    pseudo_marker = markers('pseudo'); 
    % random marker
    random_marker = markers('random');
    % target marker
    target_marker = markers('target'); 
    % before target marker
    before_target_marker = markers('before_target');
    % marker for low reaction time
    low_rt_marker = markers('low_rt');
    % marker for high reaction time
    high_rt_marker = markers('high_rt');
    % marker for HEP_target_low_RT_random_pseudo
    HEP_target_low_RT_random_pseudo_marker = markers('HEP_target_low_RT_random_pseudo');
    % marker for HEP_target_low_RT_random_pseudo
    HEP_target_high_RT_random_pseudo_marker =  markers('HEP_target_high_RT_random_pseudo');
    % marker for HEP_target_low_RT_random_pseudo
    ERP_target_low_RT_pseudo_marker =  markers('ERP_target_low_RT_random_pseudo');
    % marker for HEP_target_low_RT_random_pseudo
    ERP_target_high_RT_pseudo_marker =  markers('ERP_target_high_RT_random_pseudo');

    % Define combined markers for specific conditions
    % marker for heartbeat before target (8820)
    r_before_target_marker = [r_marker target_marker];
    % marker for target with heartbeat in the pre-stim time (8822)
    target_with_r_prestim_marker = [r_marker before_target_marker];
    % marker for pseudo heartbeat before target (48820)
    pseudo_r_before_target_marker = [pseudo_marker r_marker target_marker];
    % marker for random heartbeat before target (88820)
    random_r_before_target_marker = [random_marker r_marker target_marker];
    % marker for target with pseudo heartbeat (48822)
    target_with_pseudo_r_prestim_marker = [pseudo_marker r_marker before_target_marker];
    % marker for target with random heartbeat (88822)
    target_with_random_r_prestim_marker =[random_marker r_marker before_target_marker];

    % marker for HEP_target_low_RT (208871)
    HEP_target_low_RT_marker = [target_marker r_marker low_rt_marker];
    % marker for HEP_target_high_RT (208872)
    HEP_target_high_RT_marker = [target_marker r_marker high_rt_marker];
    % marker for ERP_target_low_RT (2071)
    ERP_target_low_RT_marker = [target_marker low_rt_marker];
    % marker for ERP_target_high_RT (2072)
    ERP_target_high_RT_marker = [target_marker high_rt_marker];

    % marker for pseudo HEP_target_low_RT (4208871)
    pseudo_HEP_target_low_RT_marker = [pseudo_marker HEP_target_low_RT_marker];
    % marker for pseudo HEP_target_high_RT (4208872)
    pseudo_HEP_target_high_RT_marker = [pseudo_marker HEP_target_high_RT_marker];
    % marker for pseudo ERP_target_low_RT (42071)
    pseudo_ERP_target_low_RT_marker = [pseudo_marker ERP_target_low_RT_marker];
    % marker for pseudo ERP_target_high_RT (42072)
    pseudo_ERP_target_high_RT_marker = [pseudo_marker ERP_target_high_RT_marker];
    
    %% Process HEP trials
    % Find indices of events matching specific markers
    idx_erp = find(ismember({EEG.event.type}, target_with_r_prestim_marker));
    idx_hep = find(ismember({EEG.event.type}, r_before_target_marker));
    RTs = [EEG.event(idx_erp).RT];
    
    % Split RTs into low and high based on median
    med_sort_val = median(RTs);
    low_RT = find(RTs<med_sort_val);
    high_RT = find(RTs>med_sort_val);

    % Backup and update event types based on RT split
    EEG_event_temp = EEG.event; 
    [EEG_event_temp(idx_hep(low_RT)).type] = deal(HEP_target_low_RT_marker);
    [EEG_event_temp(idx_hep(high_RT)).type] = deal(HEP_target_high_RT_marker);
    [EEG_event_temp(idx_erp(low_RT)).type] = deal(ERP_target_low_RT_marker);
    [EEG_event_temp(idx_erp(high_RT)).type] = deal(ERP_target_high_RT_marker);
    
    keep_indices = [idx_erp(low_RT), idx_erp(high_RT), idx_hep(low_RT), idx_hep(high_RT)];
    % Append low/high erps to original event structure, so that the original events are preserved
    EEG.event = cat(2,EEG.event,EEG_event_temp(keep_indices)); 
    
    %% Process Pseudo HEP trial
    % Find indices of events matching specific markers
    idx_erp = find(ismember({EEG.event.type}, target_with_pseudo_r_prestim_marker));
    idx_hep = find(ismember({EEG.event.type}, pseudo_r_before_target_marker));
    RTs = [EEG.event(idx_erp).RT];
    
    % Split RTs into low and high based on median
    med_sort_val = median(RTs);
    low_RT = find(RTs<med_sort_val);
    high_RT = find(RTs>med_sort_val);
    
    % Backup and update event types based on RT split
    EEG_event_temp = EEG.event;
    [EEG_event_temp(idx_hep(low_RT)).type] = deal(pseudo_HEP_target_low_RT_marker);
    [EEG_event_temp(idx_hep(high_RT)).type] = deal(pseudo_HEP_target_high_RT_marker);
    [EEG_event_temp(idx_erp(low_RT)).type] = deal(pseudo_ERP_target_low_RT_marker);
    [EEG_event_temp(idx_erp(high_RT)).type] = deal(pseudo_ERP_target_high_RT_marker);
    
    keep_indices = [idx_erp(low_RT), idx_erp(high_RT), idx_hep(low_RT), idx_hep(high_RT)];
    % Append low/high erps to original event structure, so that the original events are preserved
    EEG.event = cat(2,EEG.event,EEG_event_temp(keep_indices)); 
    
    
    %% Process random Pseudo HEP trial
    % Find indices of events matching specific markers
    idx_erp = find(ismember({EEG.event.type}, target_with_random_r_prestim_marker));
    idx_hep = find(ismember({EEG.event.type}, random_r_before_target_marker));
    RTs = [EEG.event(idx_erp).RT];
    
    % Split RTs into low and high based on median
    med_sort_val = median(RTs);
    low_RT = find(RTs < med_sort_val);
    high_RT = find(RTs > med_sort_val);
    
    % Backup and update event types based on RT split
    EEG_event_temp = EEG.event; 
    [EEG_event_temp(idx_hep(low_RT)).type] = deal(HEP_target_low_RT_random_pseudo_marker);
    [EEG_event_temp(idx_hep(high_RT)).type] = deal(HEP_target_high_RT_random_pseudo_marker);
    [EEG_event_temp(idx_erp(low_RT)).type] = deal(ERP_target_low_RT_pseudo_marker);
    [EEG_event_temp(idx_erp(high_RT)).type] = deal(ERP_target_high_RT_pseudo_marker);
    
    keep_indices = [idx_erp(low_RT), idx_erp(high_RT), idx_hep(low_RT), idx_hep(high_RT)];
    % Append low/high RT events to the original event structure
    EEG.event = cat(2, EEG.event, EEG_event_temp(keep_indices)); 
    
end

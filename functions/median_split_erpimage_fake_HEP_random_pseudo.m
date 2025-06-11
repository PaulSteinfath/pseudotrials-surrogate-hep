function [EEG, low_avg_erp, high_avg_erp, sort_val_out, idx_erp] = median_split_erpimage_fake_HEP_random_pseudo(EEG, time_window, epoch_length, filter, event_pairs, binlister_file, markers)
    % MEDIAN_SPLIT_ERPIMAGE_FAKE_HEP_RANDOM_PSEUDO - Performs median split on ERP images with pseudo-random HEP.
    %
    % Description:
    %   This function processes EEG data to perform a median split on ERP images
    %   based on a specified time window and epoch length. It applies a linear
    %   discriminant analysis (LDA) filter to the data and updates the EEG event
    %   structure with new labels for low and high amplitude target trials.
    %
    % Inputs:
    %   EEG            - EEG data structure.
    %   time_window    - Time window for HEP epoching [start, end] in milliseconds.
    %   epoch_length   - Length of the epoch in milliseconds.
    %   filter         - Linear Discriminant Analysis (LDA) filter matrix.
    %   event_pairs    - Cell array of event pairs for analysis.
    %   binlister_file - File path to the binlister file.
    %   markers        - Structure containing markers for low and high target trials.
    %
    % Outputs:
    %   EEG            - Modified EEG data structure with updated events.
    %   low_avg_erp    - Indices of trials with low average ERP.
    %   high_avg_erp   - Indices of trials with high average ERP.
    %   sort_val_out   - Sorted values used for median split.
    %   idx_erp        - Indices of ERP events in the original EEG data.
    %

    % Extract markers
    % marker for trials with low amplitude target
    low_target_marker = markers('low_target');
    % marker for trials with high amplitude target
    high_target_marker = markers('high_target');

    % Intermediate Epoching for median split
    all_events = [event_pairs{:}];
    EEGtemp = pop_selectevent(EEG, 'type', all_events, 'deleteevents', 'on');
    EEGtemp = pop_creabasiceventlist(EEGtemp, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' });
    EEGtemp = pop_binlister(EEGtemp, 'BDF', binlister_file, 'IndexEL', 1, 'SendEL2', 'EEG', 'Voutput', 'EEG', 'UpdateEEG', 'on');
    
    EEGtemp = pop_epochbin(EEGtemp, epoch_length, 'none');
    
    % Remove events that are not on latency 0
    EEGtemp = remove_non_zero_events(EEGtemp);
    
    % Loop over event pairs
    for k = 1:length(event_pairs)
        event_pair = event_pairs{k};
        
        %find new name of events after binlister renamed it
        event_task  = event_pair{1};
        even_bin_idx = find(contains({EEGtemp.event.type}, ['(',event_task,')']));
        event_bin_label = EEGtemp.event(even_bin_idx).type;
    
        % Select only data of interest
        EEGtemp_split = pop_selectevent(EEGtemp, 'type', event_bin_label, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
    
        % Check if we should apply LDA filter
        if ~isempty(filter) && size(EEGtemp_split.data, 1) >= size(filter, 2)
            % Multi-channel case with LDA
            LDA_out = [];
            for evnt = 1:size(EEGtemp_split.data, 3)
                LDA_out = [LDA_out; filter * EEGtemp_split.data(1:size(filter,2),:,evnt)];
            end
            data = LDA_out';
        else
            % Single-channel case or no LDA filter
            data = squeeze(EEGtemp_split.data(1,:,:)); % Take first channel if multiple exist
        end
    
        % Recalc time range of interest into samples
        windowsamples = find(EEGtemp_split.times >= time_window(1) & EEGtemp_split.times <= time_window(2));
    
        sort_val_out = [];
        nDim = size(data);
        for kk = 1:nDim(end)
            trial_data = data(:, kk);
            sort_val_out(kk) = mean(trial_data(windowsamples(1):windowsamples(end)));
        end
    
        % Median split into low and high
        med_sort_val = median(sort_val_out);
        low_avg_erp = find(sort_val_out < med_sort_val);
        high_avg_erp = find(sort_val_out > med_sort_val);
    
        % Take the original continuous EEG data and get indices for events of interest
        idx_erp = find(strcmp({EEG.event.type}, event_pair{1}));
        idx_hep = find(strcmp({EEG.event.type}, event_pair{2}));
    
        % Add events for the high/low epochs
        EEG_event_temp = EEG.event; % Backup of EEG.event
    
        low_event_label = strcat(event_pair{1}, low_target_marker);
        high_event_label = strcat(event_pair{1}, high_target_marker);
        low_hep_label = strcat(event_pair{2}, low_target_marker);
        high_hep_label = strcat(event_pair{2}, high_target_marker);
    
        [EEG_event_temp(idx_erp(low_avg_erp)).type] = deal(low_event_label);
        [EEG_event_temp(idx_erp(high_avg_erp)).type] = deal(high_event_label);
        [EEG_event_temp(idx_hep(low_avg_erp)).type] = deal(low_hep_label);
        [EEG_event_temp(idx_hep(high_avg_erp)).type] = deal(high_hep_label);
        
        if ~isfield(EEG.event, 'sort_val')
            [EEG.event.sort_val] = deal(NaN);
        end
    
        % Add sort_val_out to EEG event structure
        for i = 1:length(low_avg_erp)
            EEG_event_temp(idx_erp(low_avg_erp(i))).sort_val = sort_val_out(low_avg_erp(i));  % Low event sort value
            EEG_event_temp(idx_hep(low_avg_erp(i))).sort_val = sort_val_out(low_avg_erp(i));  % Low event sort value
    
        end
        for i = 1:length(high_avg_erp)
            EEG_event_temp(idx_erp(high_avg_erp(i))).sort_val = sort_val_out(high_avg_erp(i));  % High event sort value
            EEG_event_temp(idx_hep(high_avg_erp(i))).sort_val = sort_val_out(high_avg_erp(i));  % High event sort value
        end
    
        keep_indices = [idx_erp(low_avg_erp), idx_erp(high_avg_erp), idx_hep(low_avg_erp), idx_hep(high_avg_erp)];
        EEG.event = cat(2, EEG.event, EEG_event_temp(keep_indices)); % Append low/high erps to original event structure, so that the original events are preserved
    end
end

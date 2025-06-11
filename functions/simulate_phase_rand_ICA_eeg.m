function [EEG, snr_stats] = simulate_phase_rand_ICA_eeg(EEG, relationship_type, relationship_scale, sim_params, target_marker, r_marker, pseudo_marker, random_marker)
    % SIMULATE_PHASE_RAND_ICA_EEG Simulates EEG data with controlled relationships between HEP and P300
    %
    % Description:
    %   Simulates EEG data by adding P300 and HEP components with specified relationships.
    %   Allows control over signal amplitudes, timing, and relationships between components.
    %
    % Inputs:
    %   EEG               - EEGLAB EEG structure containing base data
    %   relationship_type - Type of HEP-P300 relationship ('direct', 'inverse', 'none')
    %   relationship_scale- Scale of relationship (0 to 1)
    %   sim_params        - Structure containing simulation parameters:
    %     .hep_params     - HEP parameters [mean_amp, std_amp]
    %     .p300_params    - P300 parameters [mean_amp, std_amp]
    %     .hep_timing     - HEP timing [onset, duration] in seconds
    %     .p300_timing    - P300 timing [onset, duration] in seconds
    %     .triggers       - Structure with trigger markers
    %   target_marker     - Marker for target events
    %   r_marker          - Marker for 'r' events
    %   pseudo_marker     - Marker for pseudo events
    %   random_marker     - Marker for random events
    %
    % Outputs:
    %   EEG               - Modified EEG structure with simulated components
    %   snr_stats         - Structure containing SNR statistics
    
    % Default parameters if not provided
    if ~exist('sim_params', 'var')
        sim_params = struct();
    end
    
    % Set default parameters
    default_params = struct(...
        'hep_params', [1.5, 0.8], ...     % [mean_amp, std_amp]
        'p300_params', [7, 2], ...        % [mean_amp, std_amp]
        'hep_timing', [0.2, 0.4], ...     % [onset, duration] in seconds
        'p300_timing', [0.2, 0.6], ...    % [onset, duration] in seconds
        'triggers', struct(...
            'target', target_marker, ...
            'r', r_marker, ...
            'pseudo', pseudo_marker, ...
            'random', random_marker ...
        )...
    );
    
    % Merge provided parameters with defaults
    sim_params = merge_structs(default_params, sim_params);
    
    % Extract basic EEG parameters
    npoints = EEG.pnts;
    srate = EEG.srate;
    
    % Store original data for SNR calculation
    original_EEG_data = EEG.data;
    
    % Create trigger combinations
    r_target = [sim_params.triggers('r') sim_params.triggers('target')];
    pseudo_target = [sim_params.triggers('pseudo') sim_params.triggers('target')];
    random_target = [sim_params.triggers('random') sim_params.triggers('target')];
    target_with_r_prestim_marker = [sim_params.triggers('r') sim_params.triggers('before_target')];

    % Find relevant events
    idx_p300 = find(contains({EEG.event.type}, target_with_r_prestim_marker));
    idx_hep = find(contains({EEG.event.type},r_target));
    
    % Generate amplitudes
    p300_amplitude = abs(normrnd(sim_params.p300_params(1), sim_params.p300_params(2), [1, length(idx_p300)]));
    pseudo_p300_amplitude = abs(normrnd(sim_params.p300_params(1), sim_params.p300_params(2), [1, length(idx_p300)]));
    hep_amplitude = normrnd(sim_params.hep_params(1), sim_params.hep_params(2), [1, length(idx_hep)]);
    
    % Apply relationship
    switch relationship_type
        case 'direct'
            hep_amplitude = relationship_scale * zscore(p300_amplitude) + ...
                sqrt(1 - relationship_scale^2) * zscore(hep_amplitude);
        case 'inverse'
            hep_amplitude = -relationship_scale * zscore(p300_amplitude) + ...
                sqrt(1 - relationship_scale^2) * zscore(hep_amplitude);
    end
    
    % Rescale HEP
    hep_amplitude = hep_amplitude * sim_params.hep_params(2) + sim_params.hep_params(1);
    
    % Initialize counters
    counters = struct('p300', 1, 'pseudo_p300', 1, 'pre_stim', 1);
    
    % Add components to EEG data
    for i = 1:length(EEG.event)
        event_type = EEG.event(i).type;
        
        % Add P300 for real trials
        if contains(event_type, r_target)
            EEG = add_component(EEG, i, sim_params.p300_timing, ...
                p300_amplitude(counters.p300), srate, npoints);
            counters.p300 = counters.p300 + 1;
            
        % Add HEP
        elseif strcmp(event_type, sim_params.triggers('r'))
            EEG = add_component(EEG, i, sim_params.hep_timing, ...
                hep_amplitude(counters.pre_stim), srate, npoints);
            counters.pre_stim = counters.pre_stim + 1;
            
        % Add P300 for pseudo trials
        elseif contains(event_type, pseudo_target)
            EEG = add_component(EEG, i, sim_params.p300_timing, ...
                pseudo_p300_amplitude(counters.pseudo_p300), srate, npoints);
            counters.pseudo_p300 = counters.pseudo_p300 + 1;
        end
    end
    
    % Calculate SNR statistics
    snr_stats = calculate_snr(original_EEG_data, p300_amplitude, hep_amplitude);
    
    % Store SNR in EEG structure
    EEG.snr_p300 = snr_stats.snr_p300;
    EEG.snr_HEP = snr_stats.snr_hep;
    
    EEG = eeg_checkset(EEG);
end

function run_permutations_phase_rand(EEGpath, basedir, n_subjects, n_permutations, ...
    relationship_types, relationship_scales, rm_base, baseline, sim_fs, alpha, ...
    cluster_alpha, num_randomization, channel, latency, ...
    epoch_length, stim_interval, hep_prestim_window, sim_params, settings_path, hep_epoch, lda_file, ids_file, binlister_file, markers)
% RUN_PERMUTATIONS_PHASE_RAND
% =========================================================================
% This function performs simulation and permutation analysis on preprocessed
% EEG datasets. It applies phase randomization on the designated channel, adds
% simulated triggers and aligns Heartbeat Evoked Potentials (HEP) with task events.
%
% The inputs are defined as follows (see main_preprocessing for details):
%   EEGpath             : Path to the EEG dataset files.
%   basedir             : Base directory for simulation output (sim_basedir).
%   n_subjects          : Number of subjects (sim_n_subjects).
%   n_permutations      : Number of permutation iterations (sim_n_permutations).
%   relationship_types  : Cell array containing relationship types (e.g. {'direct','inverse'}).
%   relationship_scales : Array of relationship scales (e.g. [1, 0.8, 0.6, 0.4, 0.2, 0]).
%   rm_base             : Flag indicating whether to apply baseline correction.
%   baseline            : Baseline correction interval.
%   sim_fs              : Sampling frequency for simulation.
%   alpha               : Significance level for statistics (sim_alpha).
%   cluster_alpha       : Cluster-level significance (sim_cluster_alpha).
%   num_randomization   : Number of randomizations for stats (sim_num_randomization).
%   channel             : EEG channel used for phase randomization (sim_channel).
%   latency             : Time window for statistical evaluation (sim_latency).
%   epoch_length        : Duration of epochs for ERP extraction (sim_epoch_length).
%   stim_interval       : Interval for stimulus presentation.
%   hep_prestim_window  : Time window for expected HEPs before stimulus (sim_hep_time_window).
%   sim_params          : Parameters for simulation.
%   settings_path       : Path to configuration files shared with preprocessing.
%   hep_epoch           : Time range in ms for extracting HEP events.
%   lda_file            : File with LDA filter matrices for ERP mapping.
%   ids_file            : File listing subject IDs (used to index LDA filters).
%   binlister_file      : File used for mapping events into bins (for pseudo-trial creation).
%   markers             : A containers.Map defining event markers (as in main_preprocessing).
%
% OUTPUT:
%   Statistical results and ERP plots are saved in designated directories.
% =========================================================================

%% Load configuration and set up FieldTrip defaults
load([settings_path, 'stats/layout.mat']);
load([settings_path, 'stats/neighbours.mat']);
load([settings_path, 'stats/ft_struct_250_sim_singChan.mat']);
f_struct = ft_struct;
f_struct.avg = []; f_struct.var = []; f_struct.dof = [];

ft_defaults;

%% Get EEG files and determine SLURM task ID (for clusters)
files = dir(fullfile(EEGpath, '*.set'));
arrayTaskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(arrayTaskId)
    arrayTaskId = 1;  % Default for local testing
end

%% Process each relationship type and scale
for rel_type = relationship_types
    rel_type_val = rel_type{1};
    for rel_scale = relationship_scales(arrayTaskId)
        
        % Define current output directories (only create folders that are used)
        curr_basedir   = sprintf('%s%s_scale_%g', basedir, rel_type_val, rel_scale);
        save_error     = [curr_basedir, '/Logfiles/'];
        QA_save_path   = [curr_basedir, '/QA/'];
        stats_dir      = [curr_basedir, '/stats/'];
        save_img_path  = [curr_basedir, '/stats_img/'];
        create_dirs({save_error, QA_save_path, stats_dir, save_img_path});
        
        %% Process each permutation iteration
        for perm_iter = 1:n_permutations
            if isfile([stats_dir, num2str(perm_iter), '_stat_HEP_target.mat'])
                continue;  % Skip if results already exist
            end
            
            % Initialize cell arrays for storing ERP data from each subject
            HEP_target_low              = cell(n_subjects,1);
            HEP_target_high             = cell(n_subjects,1);
            HEP_target_low_pseudo       = cell(n_subjects,1);
            HEP_target_high_pseudo      = cell(n_subjects,1);
            HEP_target_low_pseudo_2     = cell(n_subjects,1);
            HEP_target_high_pseudo_2    = cell(n_subjects,1);
            HEP_target_low_pseudo_baudry = cell(n_subjects,1);
            HEP_target_high_pseudo_baudry = cell(n_subjects,1);
            IDs = cell(n_subjects,1);
            
            %% Loop over each subject file
            for i = 1:n_subjects
                try
                    % Get subject ID from filename and create subject-specific QA dir
                    subjid = files(i).name(1:end-4);
                    create_dirs({QA_save_path}, subjid);
                    
                    % Load subject EEG data
                    EEG = pop_loadset('filename', files(i).name, 'filepath', EEGpath);
                    if size(EEG.data,1) ~= 31
                        continue;  % Process only if 31 channels exist
                    end
                    
                    %% Phase randomization (based on selected channel)
                    EEG = pop_select(EEG, 'channel', {channel});
                    EEG.data(:,1:end) = phaseran(EEG.data', 1)';
                    
                    %% Add triggers and perform baseline removal
                    EEG = add_triggers_sim(EEG, stim_interval, hep_prestim_window, markers('target_stim'), markers('r'));
                    EEG = pop_rmbase(EEG, [], []);
                    
                    %% Match HEP and task events
                    EEG = matchHEPtoTask_diffPseudo(EEG, markers('target_stim'), hep_prestim_window, ...
                        markers('r'), markers('pseudo'), markers('random'));
                    
                    %% Simulate responses with phase-randomized ICA data
                    [EEG, snr] = simulate_phase_rand_ICA_eeg(EEG, rel_type_val, rel_scale, sim_params, ...
                        markers('target'), markers('r'), markers('pseudo'), markers('random'));
                    
                    %% Adjust event markers to balance preHEP events
                    % Define combined event markers for median splitting
                    r_before_target = [markers('r'), markers('target')];
                    target_with_r_prestim  = [markers('r'), markers('before_target')];
                    pseudo_r_before_target  = [markers('pseudo'), markers('r'), markers('target')];
                    target_with_pseudo_r_prestim  = [markers('pseudo'), markers('r'), markers('before_target')];
                    random_r_before_target = [markers('random'), markers('r'), markers('target')];
                    target_with_random_r_prestim = [markers('random'), markers('r'), markers('before_target')];
                    
                    event_pairs = {
                        {target_with_r_prestim, r_before_target}, ...
                        {target_with_pseudo_r_prestim, pseudo_r_before_target}, ...
                        {target_with_random_r_prestim, random_r_before_target}
                        };
                    
                    % Calculate the number of 'r_before_target' events
                    sum_preHEP = sum(strcmp({EEG.event.type}, r_before_target));
                    % Calculate the number of 'pseudo_r_before_target' events
                    sum_pseudoPreHEP = sum(strcmp({EEG.event.type}, pseudo_r_before_target));
                    
                    % Check if the number of 'r_before_target' and 'pseudo_r_before_target' are not equal
                    if sum_preHEP ~= sum_pseudoPreHEP
                        % If there are more 'pseudo_r_before_target' events than 'r_before_target'
                        if sum_pseudoPreHEP > sum_preHEP
                            % Calculate the difference in the number of events
                            difference = sum_pseudoPreHEP - sum_preHEP;
                            % Find indices of 'pseudo_r_before_target' events
                            event_idx = find(strcmp({EEG.event.type}, pseudo_r_before_target));
                            % Randomly permute the indices
                            rem_events_idx = randperm(length(event_idx));
                            % Select the indices to be removed
                            rem_events = event_idx(rem_events_idx(1:difference));
                            % Remove the selected events from EEG
                            EEG = pop_selectevent(EEG, 'event', rem_events, 'select', 'inverse', 'deleteevents', 'on');
                            % Repeat the process for 'target_with_pseudo_r_prestim' events
                            event_idx = find(strcmp({EEG.event.type}, target_with_pseudo_r_prestim));
                            rem_events = event_idx(rem_events_idx(1:difference));
                            EEG = pop_selectevent(EEG, 'event', rem_events, 'select', 'inverse', 'deleteevents', 'on');
                            % If there are more 'r_before_target' events than 'pseudo_r_before_target'
                        elseif sum_pseudoPreHEP < sum_preHEP
                            % Calculate the difference in the number of events
                            difference = sum_preHEP - sum_pseudoPreHEP;
                            % Find indices of 'r_before_target' events
                            event_idx = find(strcmp({EEG.event.type}, r_before_target));
                            % Randomly permute the indices
                            rem_events_idx = randperm(length(event_idx));
                            % Select the indices to be removed
                            rem_events = event_idx(rem_events_idx(1:difference));
                            % Remove the selected events from EEG
                            EEG = pop_selectevent(EEG, 'event', rem_events, 'select', 'inverse', 'deleteevents', 'on');
                            % Repeat the process for 'target_with_r_prestim' events
                            event_idx = find(strcmp({EEG.event.type}, target_with_r_prestim));
                            rem_events = event_idx(rem_events_idx(1:difference));
                            EEG = pop_selectevent(EEG, 'event', rem_events, 'select', 'inverse', 'deleteevents', 'on');
                        end
                    end
                    
                    %% Load LDA filter and rescale filter weights
                    erp_t_filter = struct2array(load(lda_file));
                    ids = loadtxt(ids_file);
                    LDA_idx = find(contains(ids, subjid(1:10)));
                    filter = erp_t_filter(LDA_idx, :) * 10e-7; % Rescale from volt to microvolt
                    
                    
                    %% Split ERP images based on median (no LDA here)
                    EEG = median_split_erpimage_fake_HEP_random_pseudo(EEG, hep_epoch, epoch_length, [], event_pairs, binlister_file, markers);
                    EEG = eeg_checkset(EEG, 'eventconsistency');
                    EEG = pop_resample(EEG, sim_fs);
                    
                    %% Permute HEP timings using modified event pairs
                    % Define combined markers for specific conditions
                    r_before_low_target_marker = [markers('r'), markers('target'), markers('low_target')]; %88201
                    low_target_r_prestim_marker  = [markers('r'), markers('before_target'), markers('low_target')]; %88221
                    r_before_high_target_marker = [markers('r'), markers('target'), markers('high_target')]; %88202
                    high_target_r_prestim_marker = [markers('r'), markers('before_target'), markers('high_target')];  %88222
                    pseudo_r_before_low_target_marker = [markers('pseudo'), markers('r'), markers('target'), markers('low_target')];    %488201
                    target_with_pseudo_r_prestim_low = [markers('pseudo'), markers('r'), markers('before_target'), markers('low_target')];      %488221
                    pseudo_r_before_high_target_marker = [markers('pseudo'), markers('r'), markers('target'), markers('high_target')];     %488202
                    target_with_pseudo_r_prestim_high = [markers('pseudo'), markers('r'), markers('before_target'), markers('high_target')];   %488222
                    % Combined markers for different permutations
                    permutations_markers = markers('permutations');
                    pseudo_2_r_before_low_target_marker = [pseudo_r_before_low_target_marker permutations_markers{1}]; %4882013
                    pseudo_2_r_before_high_target_marker = [pseudo_r_before_high_target_marker permutations_markers{1}]; %4882023
                    pseudo_baudry_r_before_low_target_marker = [r_before_low_target_marker permutations_markers{1}]; %882013
                    pseudo_baudry_r_before_high_target_marker = [r_before_high_target_marker permutations_markers{1}]; %882023
                    
                    event_pairs = {
                        r_before_low_target_marker, low_target_r_prestim_marker; ...
                        r_before_high_target_marker, high_target_r_prestim_marker; ...
                        pseudo_r_before_low_target_marker, target_with_pseudo_r_prestim_low; ...
                        pseudo_r_before_high_target_marker, target_with_pseudo_r_prestim_high
                        };
                    
                    % Permutes event latencies between paired HEP and target markers
                    EEG = permute_HEP_timings(EEG,event_pairs, markers('permutations'));
                    
                    %% Create ERPs
                    % Create basic event list and apply binlister
                    EEG = pop_creabasiceventlist(EEG, 'AlphanumericCleaning', 'on', ...
                        'BoundaryNumeric', {-99}, 'BoundaryString', {'boundary'});
                    
                    EEG = pop_binlister(EEG, 'BDF', ...
                        binlister_file, ...
                        'IndexEL', 1, 'SendEL2', 'EEG', 'Voutput', 'EEG', 'UpdateEEG', 'on');
                    
                    % Create epochs
                    EEG = pop_epochbin(EEG, epoch_length, 'none');
                    EEG = remove_non_zero_events(EEG);
                    
                    % Create ERP
                    ERP = pop_averager(EEG, 'Criterion', 'good', 'SEM', 'on', 'DQ_flag', 0);
                    
                    % Baseline correction if requested
                    if rm_base
                        ERP = pop_blcerp(ERP, 'Baseline', baseline, 'Saveas', 'off');
                    end
                    
                    % Extract HEP data
                    HEP_target_low{i} = erp2fieldtrip(ERP, f_struct, r_before_low_target_marker);
                    HEP_target_high{i} = erp2fieldtrip(ERP, f_struct, r_before_high_target_marker);
                    HEP_target_low_pseudo{i} = erp2fieldtrip(ERP, f_struct, pseudo_r_before_low_target_marker);
                    HEP_target_high_pseudo{i} = erp2fieldtrip(ERP, f_struct, pseudo_r_before_high_target_marker);
                    HEP_target_low_pseudo_2{i} = erp2fieldtrip(ERP, f_struct, pseudo_2_r_before_low_target_marker);
                    HEP_target_high_pseudo_2{i} = erp2fieldtrip(ERP, f_struct, pseudo_2_r_before_high_target_marker);
                    HEP_target_low_pseudo_baudry{i} = erp2fieldtrip(ERP, f_struct, pseudo_baudry_r_before_low_target_marker);
                    HEP_target_high_pseudo_baudry{i} = erp2fieldtrip(ERP, f_struct, pseudo_baudry_r_before_high_target_marker);
                    
                    % Store subject ID
                    IDs{i} = files(i).name(1:10);
                    
                    if isfile([save_error, subjid, '_HEP_error_log.txt'])
                        delete([save_error, subjid, '_HEP_error_log.txt']);
                    end
                    
                catch ME
                    fileID = fopen([save_error, subjid, '_HEP_error_log.txt'],'w');
                    fprintf(fileID,'%6s\n',subjid);
                    fprintf(fileID,'%6s\n',getReport(ME));
                    fclose(fileID);
                end
            end
            
            % Process group results
            allsubj.ID = IDs;
            allsubj.timelock.HEP.target_low = HEP_target_low;
            allsubj.timelock.HEP.target_high = HEP_target_high;
            allsubj.timelock.HEP.target_low_pseudo = HEP_target_low_pseudo;
            allsubj.timelock.HEP.target_high_pseudo = HEP_target_high_pseudo;
            allsubj.timelock.HEP.target_low_pseudo_2 = HEP_target_low_pseudo_2;
            allsubj.timelock.HEP.target_high_pseudo_2 = HEP_target_high_pseudo_2;
            allsubj.timelock.HEP.target_low_pseudo_baudry = HEP_target_low_pseudo_baudry;
            allsubj.timelock.HEP.target_high_pseudo_baudry = HEP_target_high_pseudo_baudry;
            
            % Delete empty cells
            HEP_names = fieldnames(allsubj.timelock.HEP);
            for cond = 1:length(HEP_names)
                allsubj.timelock.HEP = rem_empty(allsubj.timelock.HEP, HEP_names(cond), IDs);
            end
            
            % Individual subject subtract fake trials
            allsubj_corr = allsubj;
            for i = 1:length(allsubj.timelock.HEP.target_high)
                allsubj_corr.timelock.HEP.target_high_pseudo{i}.avg = allsubj.timelock.HEP.target_high{i}.avg - allsubj.timelock.HEP.target_high_pseudo{i}.avg;
                allsubj_corr.timelock.HEP.target_low_pseudo{i}.avg = allsubj.timelock.HEP.target_low{i}.avg - allsubj.timelock.HEP.target_low_pseudo{i}.avg;
            end
            allsubj.timelock.HEP.target_high_correct = allsubj_corr.timelock.HEP.target_high_pseudo;
            allsubj.timelock.HEP.target_low_correct = allsubj_corr.timelock.HEP.target_low_pseudo;
            
            % Pseudotrials corrected
            allsubj_corr = allsubj;
            for i = 1:length(allsubj.timelock.HEP.target_high)
                allsubj_corr.timelock.HEP.target_high_pseudo{i}.avg = allsubj.timelock.HEP.target_high_pseudo{i}.avg - allsubj.timelock.HEP.target_high_pseudo_2{i}.avg;
                allsubj_corr.timelock.HEP.target_low_pseudo{i}.avg = allsubj.timelock.HEP.target_low_pseudo{i}.avg - allsubj.timelock.HEP.target_low_pseudo_2{i}.avg;
            end
            allsubj.timelock.HEP.target_high_pseudo_correct = allsubj_corr.timelock.HEP.target_high_pseudo;
            allsubj.timelock.HEP.target_low_pseudo_correct = allsubj_corr.timelock.HEP.target_low_pseudo;
            
            clear allsubj_corr
            
            % Move IDs out of ERP struct
            allsubj.ID = allsubj.timelock.HEP.IDs;
            allsubj.timelock.HEP = rmfield(allsubj.timelock.HEP, 'IDs');
            
            % Calculate grand averages
            cfg = [];
            cfg.channel = 'all';
            cfg.latency = latency;
            cfg.parameter = 'avg';
            cfg.keepindividual = 'no';
            
            HEP_conditions = fieldnames(allsubj.timelock.HEP);
            tempResults4 = cell(length(HEP_conditions), 1);
            parfor cond = 1:length(HEP_conditions)
                tempResults4{cond} = ft_timelockgrandaverage(cfg, allsubj.timelock.HEP.(HEP_conditions{cond}){:});
            end
            for cond = 1:length(HEP_conditions)
                allsubj.grandAvg.HEP.(HEP_conditions{cond}) = tempResults4{cond};
            end
            
            % Remove ECG
            [~, allsubj] = remove_ECG(allsubj, HEP_conditions, 'HEP');
            
            % Calculate statistics
            cfg_stats = [];
            cfg_stats.channel = channel;
            cfg_stats.latency = latency;
            cfg_stats.avgovertime = 'no';
            cfg_stats.parameter = 'avg';
            cfg_stats.neighbours = neighbours;
            cfg_stats.correctm = 'cluster';
            cfg_stats.method = 'montecarlo';
            cfg_stats.statistic = 'depsamplesT';
            cfg_stats.clusterstatistic = 'maxsum';
            cfg_stats.alpha = alpha;
            cfg_stats.clustertail = 0;
            cfg_stats.clusteralpha = cluster_alpha;
            cfg_stats.numrandomization = num_randomization;
            
            Nsub = numel(allsubj.timelock.HEP.target_low);
            cfg_stats.design(1,1:2*Nsub) = [ones(1,Nsub) 2*ones(1,Nsub)];
            cfg_stats.design(2,1:2*Nsub) = [1:Nsub 1:Nsub];
            cfg_stats.ivar = 1;
            cfg_stats.uvar = 2;
            
            % Calculate all statistical comparisons
            stat_HEP_target = ft_timelockstatistics(cfg_stats, allsubj.timelock.HEP.target_low{:}, allsubj.timelock.HEP.target_high{:});
            stat_HEP_target_pseudo = ft_timelockstatistics(cfg_stats, allsubj.timelock.HEP.target_low_pseudo{:}, allsubj.timelock.HEP.target_high_pseudo{:});
            stat_HEP_target_pseudo_baudry = ft_timelockstatistics(cfg_stats, allsubj.timelock.HEP.target_low_pseudo_baudry{:}, allsubj.timelock.HEP.target_high_pseudo_baudry{:});
            stat_HEP_target_correct = ft_timelockstatistics(cfg_stats, allsubj.timelock.HEP.target_low_correct{:}, allsubj.timelock.HEP.target_high_correct{:});
            stat_HEP_target_pseudo_correct = ft_timelockstatistics(cfg_stats, allsubj.timelock.HEP.target_low_pseudo_correct{:}, allsubj.timelock.HEP.target_high_pseudo_correct{:});
            
            % Save results
            save([stats_dir, num2str(perm_iter), '_stat_HEP_target.mat'], ...
                'stat_HEP_target', 'stat_HEP_target_correct', ...
                'stat_HEP_target_pseudo', 'stat_HEP_target_pseudo_baudry');
            
            % Generate plots
            try
                % Determine plot direction based on relationship type
                if contains(rel_type, 'direct')
                    plotDirection = 'pos';
                elseif contains(rel_type, 'inverse')
                    plotDirection = 'neg';
                else
                    plotDirection = 'pos';
                end
                
                % Plot HEP Differences
                plot_erp_var(allsubj.grandAvg.HEP.target_low, allsubj.grandAvg.HEP.target_high, ...
                    latency, {channel}, [], 'sem', stat_HEP_target, plotDirection, ...
                    {'HEP low ERP', 'HEP high ERP'}, 'northwest', [], save_img_path, ...
                    'pre_HEP_low_high_P300.png');
                
                plot_erp_var(allsubj.grandAvg.HEP.target_low_pseudo, allsubj.grandAvg.HEP.target_high_pseudo, ...
                    latency, {channel}, [], 'sem', [], plotDirection, ...
                    {'pseudo HEP low ERP', 'pseudo HEP high ERP'}, 'northeast', [], ...
                    save_img_path, 'pre_HEP_low_high_P300_pseudo.png');
                
                plot_erp_var(allsubj.grandAvg.HEP.target_low_pseudo_baudry, allsubj.grandAvg.HEP.target_high_pseudo_baudry, ...
                    latency, {channel}, [], 'sem', [], plotDirection, ...
                    {'pseudo HEP low ERP', 'pseudo HEP high ERP'}, 'northeast', [], ...
                    save_img_path, 'pre_HEP_low_high_P300_pseudo_baudry.png');
                
                plot_erp_var(allsubj.grandAvg.HEP.target_low_correct, allsubj.grandAvg.HEP.target_high_correct, ...
                    latency, {channel}, [], 'sem', stat_HEP_target_correct, plotDirection, ...
                    {'corrected \newlineHEP low ERP', 'corrected \newlineHEP high ERP'}, 'northwest', [], ...
                    save_img_path, 'tmp_pre_HEP_low_high_P300_correct.png');
                
                plot_erp_var(allsubj.grandAvg.HEP.target_low_pseudo_correct, allsubj.grandAvg.HEP.target_high_pseudo_correct, ...
                    [-1.5 1], {channel}, [], 'sem', stat_HEP_target_pseudo_correct, plotDirection, ...
                    {'low_pseudo_correct', 'high_pseudo_correct'}, 'northwest', [], ...
                    save_img_path, 'task_low_high_P300_pseudo_correct.png');
                
            catch
                warning('Failed to generate plots for permutation %d', perm_iter);
            end
            
            clear allsubj HEP_target_low HEP_target_high ...
                HEP_target_low_pseudo HEP_target_high_pseudo
            close all
        end
    end
end
end

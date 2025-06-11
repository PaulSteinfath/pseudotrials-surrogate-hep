function surrogate_statistics(data_path, save_path, save_error, binlister_file, epoch_length, lda_file, ids_file, markers, min_isi, isi, hep_prestim_window, hep_epoch, noise_epoch_tresh, baseline_correction, baseline_window, num_permutations, settings_path)
    % SURROGATE_STATISTICS - Performs permutation analysis on HEP data
    %
    % This function runs permutation tests on preprocessed heartbeat evoked potentials
    % data to analyze the relationship between HEP and behavioral/neural measures.
    %
    % Inputs:
    %   data_path        - Path to the preprocessed EEG data files (.set)
    %   save_path        - Base directory for saving outputs
    %   error_path       - Path for saving error logs
    %   qa_path          - Path for quality assurance outputs
    %   binlister_file   - File containing binlister information for ERPLab
    %   epoch_length     - Length of epoch for processing [start_ms, end_ms]
    %   lda_file         - File containing LDA matrices
    %   ids_file         - File containing subject IDs
    %   markers          - Structure containing event markers
    %   min_isi          - Minimum inter-stimulus interval allowed (seconds)
    %   isi              - Target inter-stimulus interval (seconds)
    %   hep_prestim_window - Pre-stimulus window for HEP analysis [start_s, end_s]
    %   hep_epoch        - HEP epoch length [start_ms, end_ms]
    %   noise_epoch_tresh - Noise threshold for epoch rejection (µV)
    %   baseline_correction - Whether to perform baseline correction (1=yes, 0=no)
    %   baseline_window  - Window for baseline correction [start_ms, end_ms]
    %   num_permutations - Number of permutations to perform
    %
    % Notes:
    %   - The function uses parallel processing for faster runtime during permutation
    %     testing and subject-level data processing.
    %   - Requires the FieldTrip toolbox to be installed and properly set up.
    
    %% Add required paths
    addpath(fileparts(mfilename('fullpath'))) % Add functions directory
    
    % Check for fieldtrip
    if ~exist('ft_defaults', 'file')
        error('FieldTrip not found. Please make sure it is installed and on the MATLAB path.');
    end
    ft_defaults
   
    %% Adjust output directories based on baseline correction
    if baseline_correction == 1
        base_folder = 'baseline';
    else
        base_folder = 'no_baseline';
    end
    save_path = fullfile(save_path, base_folder);
    
    %% Create output directories with prefix
    stats_dir = fullfile(save_path,'stats');
    save_img_path = fullfile(save_path,  'stats_img');
    
    create_dirs({stats_dir, save_img_path});
         
    %% Prepare FieldTrip structures
 
    % Load FieldTrip structures
    load(fullfile(settings_path, 'stats/layout.mat'));
    load(fullfile(settings_path, 'stats/neighbours.mat')); 
    load(fullfile(settings_path, 'stats/f_struct_250_3s_32chan.mat'));
    f_struct.avg = [];  
    f_struct.var = []; 
    f_struct.dof = [];
    
    %% Get files
    files = dir(fullfile(data_path, '*.set'));
    
    %% Minimum number of events per subject
    % Exclude subject if it has fewer than min_event_nb events of interest
    min_event_nb = 20;

    %% Load LDA matrices
    erp_t_filter = struct2array(load(lda_file));
    ids = loadtxt(ids_file);
    
    %% Run permutations
    for perm_iter = 1:num_permutations
        % Randomize seed for permutations
        rng('shuffle');
        
        % Skip if already done
        if isfile(fullfile(stats_dir, [num2str(perm_iter), '_stat_HEP_target_RT.mat'])) || ...
           isfile(fullfile(stats_dir, [num2str(perm_iter), '_stat_HEP_target.mat']))
            continue
        end
        
        % Initialize arrays for storing results
        HEP_target_low = cell(length(files), 1);
        HEP_target_high = cell(length(files), 1);
        HEP_target_low_RT = cell(length(files), 1);
        HEP_target_high_RT = cell(length(files), 1);
        IDs = cell(length(files), 1);
        
        % Process each subject - 20 for testing, otherwise use all
        % parfor i = 1:20 % if testing
        parfor i = 1:length(files) 
            try
                % Get subject ID
                subjid = files(i).name(1:end-4);
                                
                % Load data
                EEG = pop_loadset('filename', files(i).name, 'filepath', data_path);
                
                % De-mean the signal
                EEG = pop_rmbase(EEG, [], []);
                
                % Check if task events are present
                EEG = check_events(EEG, markers('target_stim'), markers('response'), min_isi, save_error, subjid);
                
                % Detect R-Peaks and add them to EEG
                [EEG, ~] = insert_R_peaks(EEG, EEG.ECG.data(1,:), 1, markers('r'));
                
                % Add reaction time
                EEG = calcRT(EEG, markers('target_stim'), markers('response'), isi);
                
                % Match HEP to task
                EEG = matchHEPtoTask_diffPseudo(EEG, markers('target_stim'), hep_prestim_window, markers('r'), markers('pseudo'), markers('random'));
                
                % Split by RT
                EEG = split_by_RT_pseudoHEP(EEG, markers);
                EEG = split_by_RT_target(EEG, markers);
                
                % Check events consistency
                EEG = eeg_checkset(EEG, 'eventconsistency');
                
                % Split by activity
                % Match LDA id's with EEG id's
                LDA_idx = find(contains(ids, subjid(1:10)));
                filter = erp_t_filter(LDA_idx, :) * 10e-7; % Rescale filter weights
                
               
                % Define combined markers for specific conditions
                % Basic markers (from stats_sort_random_pseudo.m)
                r_before_target_marker = [markers('r') markers('target')];
                target_with_r_prestim_marker = [markers('r') markers('before_target')];
                pseudo_r_before_target_marker = [markers('pseudo') markers('r') markers('target')];
                target_with_pseudo_r_prestim_marker = [markers('pseudo') markers('r') markers('before_target')];
                random_r_before_target_marker = [markers('random') markers('r') markers('target')];
                target_with_random_r_prestim_marker = [markers('random') markers('r') markers('before_target')];

                % Additional markers for RT and target amplitude conditions
                r_before_low_target_marker = [r_before_target_marker markers('low_target')];    % 88201
                r_before_high_target_marker = [r_before_target_marker markers('high_target')];  % 88202
                target_low_r_prestim_marker = [target_with_r_prestim_marker markers('low_target')];    % 88221
                target_high_r_prestim_marker = [target_with_r_prestim_marker markers('high_target')];  % 88222
                
                % RT-specific markers
                HEP_target_low_RT_marker = [markers('target') markers('r') markers('low_rt')];   % 208871
                HEP_target_high_RT_marker = [markers('target') markers('r') markers('high_rt')];  % 208872
                ERP_target_low_RT_marker = [markers('target') markers('low_rt')];    % 2071
                ERP_target_high_RT_marker = [markers('target') markers('high_rt')];   % 2072

     
                event_pairs = {
                    {target_with_r_prestim_marker, r_before_target_marker}, ...
                    {target_with_pseudo_r_prestim_marker, pseudo_r_before_target_marker}, ...
                    {target_with_random_r_prestim_marker, random_r_before_target_marker},
                };
                
                % Using LDA
                EEG = median_split_erpimage_fake_HEP_random_pseudo(EEG, hep_epoch, epoch_length, filter, event_pairs, binlister_file, markers);
                
                % Check event consistency
                EEG = eeg_checkset(EEG, 'eventconsistency');
                
                event_pairs = {
                    r_before_high_target_marker, target_high_r_prestim_marker;  % 88202 - 88222
                    r_before_low_target_marker, target_low_r_prestim_marker;    % 88201 - 88221
                    HEP_target_low_RT_marker, ERP_target_low_RT_marker;        % 208871 - 2071
                    HEP_target_high_RT_marker, ERP_target_high_RT_marker;      % 208872 - 2072
                };

                % Permute timings per condition
                EEG = permute_HEP_timings(EEG, event_pairs, markers('permutations'));
                
                % Prepare data for epoching in ERPLab format
                EEG = pop_creabasiceventlist(EEG, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', {-99}, 'BoundaryString', {'boundary'});
                EEG = pop_binlister(EEG, 'BDF', binlister_file, 'IndexEL', 1, 'SendEL2', 'EEG', 'Voutput', 'EEG', 'UpdateEEG', 'on');
                
                % Epoch without baseline correction
                EEG = pop_epochbin(EEG, [epoch_length(1) epoch_length(2)], 'none');
                
                % Clean epochs
                EEG = remove_non_zero_events(EEG);
                
                % Remove noisy epochs
                EEG = pop_artmwppth(EEG, 'Channel', 1:31, 'Flag', 1, 'Threshold', noise_epoch_tresh, 'Twindow', [epoch_length(1) epoch_length(2)], 'Windowsize', 200, 'Windowstep', 100);
                EEG = pop_rejepoch(EEG, find(EEG.reject.rejmanual), 0);
                
                % Create average ERP set
                ERP = pop_averager(EEG, 'Criterion', 'good', 'DQ_flag', 1, 'DSindex', 1, 'ExcludeBoundary', 'on', 'SEM', 'on');
         
                % get the list of events
                eventList = {ERP.EVENTLIST.bdf.expression};
        
                % how many random pseudotrials for this subject?
                pre_target_random_pseudo_idx = find(strcmp(eventList, ['.{' r_before_target_marker '}']));
                % pre_target_random_pseudo_idx = find(strcmp(eventList, '.{88820}'));
                % If fewer than min_event_nb ERPs, subject is excluded
                if ERP.EVENTLIST.trialsperbin(pre_target_random_pseudo_idx)  <= min_event_nb
                    continue
                end
                        
                % Apply baseline correction if needed
                if baseline_correction == 1
                    ERP = pop_blcerp(ERP, 'Baseline', baseline_window, 'Saveas', 'off');
                end
                
                perm_markers = markers('permutations');

                % Additional markers for permuted trials
                perm_target_low_r_marker = [r_before_low_target_marker perm_markers(1)];      % 882013 
                perm_target_high_r_marker = [r_before_high_target_marker perm_markers(1)];    % 882023
                perm_target_low_RT_marker = [HEP_target_low_RT_marker perm_markers(1)];       % 2088713
                perm_target_high_RT_marker = [HEP_target_high_RT_marker perm_markers(1)];     % 2088723

                % Convert to fieldtrip format using the combined markers
                HEP_target_low{i} = erp2fieldtrip(ERP, f_struct, [perm_target_low_r_marker{:}]);
                HEP_target_high{i} = erp2fieldtrip(ERP, f_struct, [perm_target_high_r_marker{:}]); 
                HEP_target_low_RT{i} = erp2fieldtrip(ERP, f_struct, [perm_target_low_RT_marker{:}]);
                HEP_target_high_RT{i} = erp2fieldtrip(ERP, f_struct, [perm_target_high_RT_marker{:}]);
                
                % Store subject ID
                IDs{i} = files(i).name(1:10);
                
                % Clear error logs if successful
                error_log_file = fullfile(save_error, [subjid, '_HEP_error_log.txt']);
                if isfile(error_log_file)
                    delete(error_log_file);
                end
                
            catch ME
                fileID = fopen(fullfile(save_error, [subjid, '_HEP_error_log.txt']), 'w');
                fprintf(fileID, '%6s\n', subjid);
                fprintf(fileID, '%6s\n', getReport(ME));
                fclose(fileID);
            end
        end
        
        % Store results
        allsubj.ID = IDs;
        allsubj.timelock.HEP.target_low = HEP_target_low;
        allsubj.timelock.HEP.target_high = HEP_target_high;
        allsubj.timelock.HEP.target_low_RT = HEP_target_low_RT;
        allsubj.timelock.HEP.target_high_RT = HEP_target_high_RT;
        
        % Delete empty cells and bad epochs
        HEP_names = fieldnames(allsubj.timelock.HEP);
        for cond =1:length(HEP_names)
            allsubj.timelock.HEP = rem_empty(allsubj.timelock.HEP, HEP_names(cond), IDs);
        end
        
        % Move IDs out of ERP struct
        allsubj.ID = allsubj.timelock.HEP.IDs;
        allsubj.timelock.HEP = rmfield(allsubj.timelock.HEP, 'IDs');
        
        % Backup
        task_allsubj = allsubj;
        
        % Grand average
        cfg = [];
        cfg.channel = 'all';
        cfg.parameter = 'avg';
        cfg.keepindividual = 'no';
        cfg.latency = [-0.2 0.6];
        
        HEP_conditions = fieldnames(task_allsubj.timelock.HEP);
        tempResults = cell(length(HEP_conditions), 1);
        
        parfor cond = 1:length(HEP_conditions)
            tempResults{cond} = ft_timelockgrandaverage(cfg, task_allsubj.timelock.HEP.(HEP_conditions{cond}){:});
        end
        
        for cond = 1:length(HEP_conditions)
            task_allsubj.grandAvg.HEP.(HEP_conditions{cond}) = tempResults{cond};
        end
        
        clear tempResults;
        
        % Remove ECG
        HEP_conditions = fieldnames(task_allsubj.timelock.HEP);
        [ECG_rest_allsubj_HEP, task_allsubj] = remove_ECG(task_allsubj, HEP_conditions, 'HEP');
        
        % Permutation test settings
        cfg_stats = [];
        cfg_stats.channel = ['all'];
        cfg_stats.latency = [-0.2 0.6];
        cfg_stats.avgovertime = 'no';
        cfg_stats.parameter = 'avg';
        cfg_stats.neighbours = neighbours;
        cfg_stats.correctm = 'cluster';
        cfg_stats.method = 'montecarlo';
        cfg_stats.statistic = 'depsamplesT';
        cfg_stats.clusterstatistic = 'maxsum';
        cfg_stats.alpha = 0.05;
        cfg_stats.clustertail = 0;
        cfg_stats.clusteralpha = 0.05;
        cfg_stats.numrandomization = 1000; % reduce for testing: 100
        
        Nsub = numel(task_allsubj.timelock.HEP.target_low);
        cfg_stats.design(1,1:2*Nsub) = [ones(1,Nsub) 2*ones(1,Nsub)];
        cfg_stats.design(2,1:2*Nsub) = [1:Nsub 1:Nsub];
        cfg_stats.ivar = 1;
        cfg_stats.uvar = 2;
        cfg_stats.minnbchan = 2;
        
        % Run the statistics
        stat_HEP_target = ft_timelockstatistics(cfg_stats, task_allsubj.timelock.HEP.target_low{:}, task_allsubj.timelock.HEP.target_high{:});
        stat_HEP_target_RT = ft_timelockstatistics(cfg_stats, task_allsubj.timelock.HEP.target_low_RT{:}, task_allsubj.timelock.HEP.target_high_RT{:});
        
        % Save results
        save(fullfile(stats_dir, [num2str(perm_iter), '_stat_HEP_target.mat']), 'stat_HEP_target');
        save(fullfile(stats_dir, [num2str(perm_iter), '_stat_HEP_target_RT.mat']), 'stat_HEP_target_RT');
        
        % Generate visualization plots
        cfg = [];
        cfg.operation = 'subtract';
        cfg.parameter = 'avg';
        
        % Difference: low - high ERP
        diff_HEP_target = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low, task_allsubj.grandAvg.HEP.target_high);
        diff_HEP_target_RT = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low_RT, task_allsubj.grandAvg.HEP.target_high_RT);
        
        % Plot topographies
        cfg_topo = [];
        cfg_topo.channel = 1:31; % Exclude ECG channel
        cfg_topo.highlight = 'on';
        cfg_topo.comment = 'no';
        cfg_topo.selchan = stat_HEP_target.label;
        cfg_topo.layout = layout;
        cfg_topo.interactive = 'no';
        
        % Truncate labels if needed
        if length(diff_HEP_target.label) > 31
            diff_HEP_target.label = diff_HEP_target.label(1:31);
        end
        
        % First plot
        cfg_topo = plot_sig_topo(stat_HEP_target, cfg_topo, diff_HEP_target, diff_HEP_target.time, 1, 'pos', 0);
        print('-dsvg', fullfile(save_img_path, [num2str(perm_iter), 'P300_HEP_topo.svg']));
        
        % Second plot (RT)
        if length(diff_HEP_target_RT.label) > 31
            diff_HEP_target_RT.label = diff_HEP_target_RT.label(1:31);
        end
        cfg_topo_RT = plot_sig_topo(stat_HEP_target_RT, cfg_topo, diff_HEP_target_RT, diff_HEP_target_RT.time, 1, 'neg', 0);
        print('-dsvg', fullfile(save_img_path, [num2str(perm_iter), 'RT_HEP_topo.svg']));
        
        % Plot ERPs using highlighted channels if significant clusters exist, otherwise use Pz
        if isfield(stat_HEP_target, 'posclusters') && ~isempty(stat_HEP_target.posclusters) && ...
           any([stat_HEP_target.posclusters.prob] <= 0.025)
            plot_channels = cfg_topo.highlightchannel;
        else
            plot_channels = {'Pz'};
        end
        
        plot_erp_var(task_allsubj.grandAvg.HEP.target_low, task_allsubj.grandAvg.HEP.target_high, ...
                    [-0.2 0.6], plot_channels, [], 'sem', stat_HEP_target, 'pos', ...
                    {'HEP low ERP', 'HEP high ERP'}, 'northwest', [], save_img_path, ...
                    [num2str(perm_iter), '_pre_HEP_low_high_P300.png']);
        
        % Same for RT ERPs
        if isfield(stat_HEP_target_RT, 'posclusters') && ~isempty(stat_HEP_target_RT.posclusters) && ...
           any([stat_HEP_target_RT.posclusters.prob] <= 0.025)
            plot_channels_RT = cfg_topo_RT.highlightchannel;
        else
            plot_channels_RT = {'Pz'};
        end
        
        plot_erp_var(task_allsubj.grandAvg.HEP.target_low_RT, task_allsubj.grandAvg.HEP.target_high_RT, ...
                    [-0.2 0.6], plot_channels_RT, [], 'sem', stat_HEP_target_RT, 'neg', ...
                    {'HEP fast RT', 'HEP slow RT'}, 'northwest', [],save_img_path, ...
                    [num2str(perm_iter), '_pre_HEP_low_high_RT.png']);
                
        close all;
        clear task_allsubj allsubj;
    end
end

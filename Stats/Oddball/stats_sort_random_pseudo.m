% Script for the following statistical analyses of HEPs and ERPs
% for Task and Rest data
% Includes:
% - optional baseline correction for HEP
% - selectiom of relevant events
% - cluster-based permutation t-tests (two-tailed) using the
% FieldTrip toolbox to compare HEPs between conditions
% The script also applies statistical masks to the data,
% saves the results, creates contrasts, and generates topographical plots for significant clusters.

%% Task or Rest data
% Perform statistics on task or rest data
% Important: Be sure that the variable task is set in the main_preprocessing script
% task = 0; % 0: rest, 1: task
task = 1; % 0: rest, 1: task

%% Paths
if task == 1
    data_path = task_epochs_path;
%    data_path = 'D:\Code_review_Paul\example_data_epochs\Task/';
    save_path = task_output_path; 
    error_path = task_error_path;
else
    data_path = rest_epochs_path;
%    data_path = 'D:\Code_review_Paul\example_data_epochs\Rest/';
    save_path = rest_output_path; 
    error_path = rest_error_path;
end

%% Minimum number of events per subject
% Exclude subject if it has fewer than min_event_nb events of interest
min_event_nb = 20;

%% Markers
% Define combined markers for specific conditions
% marker for heartbeat before target (8820)
r_before_target_marker = [markers('r') markers('target')];
% marker for random heartbeat before target (88820)
random_r_before_target_marker = [markers('random') markers('r') markers('target')];
% marker for target with random heartbeat (88822)
target_with_random_r_prestim_marker =[markers('random') markers('r') markers('before_target')];

% marker for low target with random heartbeat (888221)
low_target_with_random_r_prestim_marker = [target_with_random_r_prestim_marker markers('low_target')];
% marker for high target with random heartbeat (888222)
high_target_with_random_r_prestim_marker = [target_with_random_r_prestim_marker markers('high_target')];
% marker for heartbeat before low target (88201)
r_before_low_target_marker = [r_before_target_marker markers('low_target')];
% marker for heartbeat before high target (88202)
r_before_high_target_marker = [r_before_target_marker markers('high_target')];
% marker for random heartbeat before low target (888201)
random_r_before_low_target_marker = [random_r_before_target_marker markers('low_target')];
% marker for random heartbeat before high target (888202)
random_r_before_high_target_marker = [random_r_before_target_marker markers('high_target')];

% marker for HEP_target_low_RT (208871)
HEP_target_low_RT_marker = [target_marker r_marker low_rt_marker];
% marker for HEP_target_high_RT (208872)
HEP_target_high_RT_marker = [target_marker r_marker high_rt_marker];
% marker for ERP_target_low_RT (2071)
ERP_target_low_RT_marker = [target_marker low_rt_marker];
% marker for ERP_target_high_RT (2072)
ERP_target_high_RT_marker = [target_marker high_rt_marker];

% marker for HEP_target_low_RT_random_pseudo
HEP_target_low_RT_random_pseudo_marker = markers('HEP_target_low_RT_random_pseudo');
% marker for HEP_target_low_RT_random_pseudo
HEP_target_high_RT_random_pseudo_marker =  markers('HEP_target_high_RT_random_pseudo');
% marker for HEP_target_low_RT_random_pseudo
ERP_target_low_RT_pseudo_marker =  markers('ERP_target_low_RT_random_pseudo');
% marker for HEP_target_low_RT_random_pseudo
ERP_target_high_RT_pseudo_marker =  markers('ERP_target_high_RT_random_pseudo');


%% files to process
files = dir([data_path, '*.erp']);

%% load included subjects
incl_list = readtable([settings_path, 'inclusion_exclusion/included_subj.xlsx']);

%% remove additional subjects from list
load([settings_path, 'inclusion_exclusion/remove_subjects.mat'])
% remove subjects with bad ECG
ECG_remove = ismember(incl_list.incl_list, bad_ECG_ids);
% remove subjects with bad HEP
HEP_remove = ismember(incl_list.incl_list, bad_HEP_ids);
rows_to_remove = ECG_remove | HEP_remove;
incl_list(rows_to_remove, :) = [];

%% save path
if rm_base == 1
    savepath = [save_path, 'baseline/'];
else 
    savepath = [save_path, 'no_baseline/'];
end
create_dirs({savepath});

%% load some prerequisites 
% EEG layout and electrodes location
load([settings_path, 'stats/layout.mat'])
load([settings_path, 'stats/neighbours.mat'])
load([settings_path, 'stats/f_struct_250_3s_32chan.mat'])
f_struct.avg = [];  f_struct.var = []; f_struct.dof = [];

%% Initialize cell arrays
numFiles = length(files);
ERP_target_low = cell(1, numFiles);
ERP_target_high = cell(1, numFiles);
HEP_target_low = cell(1, numFiles);
HEP_target_high = cell(1, numFiles);
HEP_target_low_random_pseudo = cell(1, numFiles);
HEP_target_high_random_pseudo = cell(1, numFiles);
HEP_target_low_RT = cell(1, numFiles);
HEP_target_high_RT = cell(1, numFiles);
ERP_target_low_RT = cell(1, numFiles);
ERP_target_high_RT = cell(1, numFiles);
HEP_target_low_RT_random_pseudo = cell(1, numFiles);
HEP_target_high_RT_random_pseudo = cell(1, numFiles);
ERP_target_low_RT_pseudo = cell(1, numFiles);
ERP_target_high_RT_pseudo = cell(1, numFiles);
HEP_target_random_pseudo = cell(1, numFiles);
IDs = cell(1, numFiles);

% loop over subjects
for subj = 1:length(files)
% parfor subj = 1:length(files)
    try

        % get subjid
        subjid = files(subj).name(1:10);

        % check if on incl_list. if not, skip it
        if find(contains(table2cell(incl_list), subjid)) == 0
            continue
        end

        % load data
        ERP = pop_loaderp('filename',[files(subj).name(1:end-4), '.erp'],'filepath',data_path,'overwrite','off','Warning','on');

        % get the list of events
        eventList = {ERP.EVENTLIST.bdf.expression};

        % how many random pseudotrials for this subject?
        pre_target_random_pseudo_idx = find(strcmp(eventList, ['.{' random_r_before_target_marker '}']));
        % pre_target_random_pseudo_idx = find(strcmp(eventList, '.{88820}'));
        % If fewer than min_event_nb ERPs, subject is excluded
        if ERP.EVENTLIST.trialsperbin(pre_target_random_pseudo_idx)  <= min_event_nb
            continue
        end

        % how many HEPs for this subject?
        pre_target_heartbeat_idx = find(strcmp(eventList, ['.{' r_before_target_marker '}']));
        % pre_target_heartbeat_idx = find(strcmp(eventList, '.{8820}'));
        % if fewer than min_event_nb pre-stimulus HEPs, subject is excluded
        if ERP.EVENTLIST.trialsperbin(pre_target_heartbeat_idx)  <= min_event_nb
            continue
        end


        %% baseline correction
        if rm_base == 1
            ERP = pop_blcerp(ERP, 'Baseline', baseline, 'Saveas', 'off' );
        end

        %% High Low ERP
        ERP_target_low{subj} =   erp2fieldtrip(ERP, f_struct, low_target_with_random_r_prestim_marker);
        ERP_target_high{subj} =  erp2fieldtrip(ERP, f_struct, high_target_with_random_r_prestim_marker);
        % ERP_target_low{subj} =   erp2fieldtrip(ERP, f_struct, '88221');
        % ERP_target_high{subj} =  erp2fieldtrip(ERP, f_struct, '88222');

        %% High Low HEP
        HEP_target_low{subj} =  erp2fieldtrip(ERP, f_struct, r_before_low_target_marker);
        HEP_target_high{subj} =  erp2fieldtrip(ERP, f_struct, r_before_high_target_marker);
        % HEP_target_low{subj} =  erp2fieldtrip(ERP, f_struct, '88201');
        % HEP_target_high{subj} =  erp2fieldtrip(ERP, f_struct, '88202');

        %% High Low random pseudo HEP
        HEP_target_low_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, random_r_before_low_target_marker);
        HEP_target_high_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, random_r_before_high_target_marker);
        % HEP_target_low_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, '888201');
        % HEP_target_high_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, '888202');

        %% High Low RT
        HEP_target_low_RT{subj} =  erp2fieldtrip(ERP, f_struct, HEP_target_low_RT_marker);
        HEP_target_high_RT{subj} =  erp2fieldtrip(ERP, f_struct, HEP_target_high_RT_marker);
        ERP_target_low_RT{subj} =  erp2fieldtrip(ERP, f_struct, ERP_target_low_RT_marker);
        ERP_target_high_RT{subj} =  erp2fieldtrip(ERP, f_struct, ERP_target_high_RT_marker);
        % HEP_target_low_RT{subj} =  erp2fieldtrip(ERP, f_struct, '208871');
        % HEP_target_high_RT{subj} =  erp2fieldtrip(ERP, f_struct, '208872');
        % ERP_target_low_RT{subj} =  erp2fieldtrip(ERP, f_struct, '2071';
        % ERP_target_high_RT{subj} =  erp2fieldtrip(ERP, f_struct, '2072');

        %% High Low random Pseudo RT
        HEP_target_low_RT_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, HEP_target_low_RT_random_pseudo_marker);
        HEP_target_high_RT_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, HEP_target_high_RT_random_pseudo_marker);
        ERP_target_low_RT_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, ERP_target_low_RT_pseudo_marker);
        ERP_target_high_RT_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, ERP_target_high_RT_pseudo_marker);
        % HEP_target_low_RT_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, '15');
        % HEP_target_high_RT_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, '16');
        % ERP_target_low_RT_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, '17');
        % ERP_target_high_RT_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, '18');

        %% random pseudo HEP
        HEP_target_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, random_r_before_target_marker);
        % HEP_target_random_pseudo{subj} =  erp2fieldtrip(ERP, f_struct, '88820');

        %% Store SubjID
        IDs{subj} = files(subj).name(1:10);

    catch ME
        % if error, save it
        fileID = fopen([error_path, subjid, '_stats_error_log.txt'],'w');
        fprintf(fileID,'%6s\n',subjid);
        fprintf(fileID,'%6s\n',getReport(ME));
        fclose(fileID);
    end
end

%% Put in struct for fieldtrip
task_allsubj = {};
task_allsubj.ID = IDs;

task_allsubj.timelock.ERP.target_low = ERP_target_low;
task_allsubj.timelock.ERP.target_high = ERP_target_high;

task_allsubj.timelock.HEP.target_low = HEP_target_low;
task_allsubj.timelock.HEP.target_high = HEP_target_high;

task_allsubj.timelock.ERP.target_low_RT = ERP_target_low_RT;
task_allsubj.timelock.ERP.target_high_RT = ERP_target_high_RT;
task_allsubj.timelock.HEP.target_low_RT = HEP_target_low_RT;
task_allsubj.timelock.HEP.target_high_RT = HEP_target_high_RT;

task_allsubj.timelock.HEP.target_low_random_pseudo = HEP_target_low_random_pseudo;
task_allsubj.timelock.HEP.target_high_random_pseudo = HEP_target_high_random_pseudo;

task_allsubj.timelock.HEP.target_low_RT_random_pseudo = HEP_target_low_RT_random_pseudo;
task_allsubj.timelock.HEP.target_high_RT_random_pseudo = HEP_target_high_RT_random_pseudo;

task_allsubj.timelock.HEP.target_random_pseudo = HEP_target_random_pseudo;

%% delete empty cells 
ERP_names = fieldnames(task_allsubj.timelock.ERP);
for cond =1:length(ERP_names)
    task_allsubj.timelock.ERP = rem_empty(task_allsubj.timelock.ERP, ERP_names(cond), IDs);
end

HEP_names = fieldnames(task_allsubj.timelock.HEP);
for cond =1:length(HEP_names)
    task_allsubj.timelock.HEP = rem_empty(task_allsubj.timelock.HEP, HEP_names(cond), IDs);
end

%% Move IDs out of ERP and HEP struct into ID
task_allsubj.ID  = task_allsubj.timelock.HEP.IDs;
task_allsubj.timelock.HEP = rmfield(task_allsubj.timelock.HEP,'IDs');
task_allsubj.timelock.ERP = rmfield(task_allsubj.timelock.ERP,'IDs');

%% clear some memory
clear ERP_standard ERP_target ERP_novelty ERP_target_low ERP_target_high ERP_standard_low ERP_standard_high HEP_target_low HEP_target_high HEP_standard_low ...
    HEP_standard_high ERP_target_low_RT ERP_target_high_RT HEP_target_low_RT HEP_target_low_RT HEP_target_high_RT

%% Individual subject subtract fake trials
task_allsubj_corr = task_allsubj;

% task corrected random pseudo
for i = 1:length(task_allsubj.timelock.HEP.target_random_pseudo)
    % Copy full structure and only modify the avg field 
    task_allsubj_corr.timelock.HEP.target_high_pseudo{i} = task_allsubj.timelock.HEP.target_high{i};
    task_allsubj_corr.timelock.HEP.target_low_pseudo{i} = task_allsubj.timelock.HEP.target_low{i};
    task_allsubj_corr.timelock.HEP.target_low_RT_pseudo{i} = task_allsubj.timelock.HEP.target_low_RT{i};
    task_allsubj_corr.timelock.HEP.target_high_RT_pseudo{i} = task_allsubj.timelock.HEP.target_high_RT{i};
    
    % Subtract averages
    task_allsubj_corr.timelock.HEP.target_high_pseudo{i}.avg = task_allsubj.timelock.HEP.target_high{i}.avg - task_allsubj.timelock.HEP.target_high_random_pseudo{i}.avg;
    task_allsubj_corr.timelock.HEP.target_low_pseudo{i}.avg = task_allsubj.timelock.HEP.target_low{i}.avg - task_allsubj.timelock.HEP.target_low_random_pseudo{i}.avg;
    task_allsubj_corr.timelock.HEP.target_low_RT_pseudo{i}.avg = task_allsubj.timelock.HEP.target_low_RT{i}.avg - task_allsubj.timelock.HEP.target_low_RT_random_pseudo{i}.avg;
    task_allsubj_corr.timelock.HEP.target_high_RT_pseudo{i}.avg = task_allsubj.timelock.HEP.target_high_RT{i}.avg - task_allsubj.timelock.HEP.target_high_RT_random_pseudo{i}.avg;
end

% Store corrected random data
task_allsubj.timelock.HEP.target_high_correct_random = task_allsubj_corr.timelock.HEP.target_high_pseudo;
task_allsubj.timelock.HEP.target_low_correct_random = task_allsubj_corr.timelock.HEP.target_low_pseudo;
task_allsubj.timelock.HEP.target_low_RT_correct_random = task_allsubj_corr.timelock.HEP.target_low_RT_pseudo;
task_allsubj.timelock.HEP.target_high_RT_correct_random = task_allsubj_corr.timelock.HEP.target_high_RT_pseudo;

clear task_allsubj_corr

%% perform grand averaging
% settings for grand averaging
cfg_ga = [];
cfg_ga.channel = 'all';
cfg_ga.latency = [-0.2 0.6];
cfg_ga.parameter = 'avg';
cfg_ga.keepindividual = 'no';

% for ERPs
ERP_conditions = fieldnames(task_allsubj.timelock.ERP);
tempResults1 = cell(length(ERP_conditions), 1);
parfor cond = 1:length(ERP_conditions)
    tempResults1{cond} = ft_timelockgrandaverage(cfg_ga, task_allsubj.timelock.ERP.(ERP_conditions{cond}){1,:});
end
for cond = 1:length(ERP_conditions)
    task_allsubj.grandAvg.ERP.(ERP_conditions{cond}) = tempResults1{cond};
end

% for HEPs
HEP_conditions = fieldnames(task_allsubj.timelock.HEP);
tempResults = cell(length(HEP_conditions), 1);
parfor cond  = 1:length(HEP_conditions)
    tempResults{cond} = ft_timelockgrandaverage(cfg_ga, task_allsubj.timelock.HEP.(HEP_conditions{cond}){1,:});
end
for cond = 1:length(HEP_conditions)
    task_allsubj.grandAvg.HEP.(HEP_conditions{cond}) = tempResults{cond};
end

%% Take out the ECG to own struct
ERP_conditions = fieldnames(task_allsubj.timelock.ERP);
[ECG_task_allsubj, task_allsubj] = remove_ECG(task_allsubj, ERP_conditions, 'ERP');

HEP_conditions = fieldnames(task_allsubj.timelock.HEP);
[ECG_task_allsubj_HEP, task_allsubj] = remove_ECG(task_allsubj, HEP_conditions, 'HEP');

%% permutation settings
cfg_stats = [];
cfg_stats.channel     = ['all'];
cfg_stats.latency     = [];
cfg_stats.avgovertime = 'no';
cfg_stats.parameter   = 'avg';
cfg_stats.neighbours       = neighbours;  
cfg_stats.correctm         = 'cluster';
cfg_stats.method      = 'montecarlo';
cfg_stats.statistic   = 'depsamplesT';
cfg_stats.clusterstatistic = 'maxsum';

%Keeping alpha at 0.05 for one-sided t-tests 
cfg_stats.alpha       = 0.05;

cfg_stats.correcttail = 'alpha'; % no need at this stage
cfg_stats.clustertail      = 0;

%take a initial significance threshold of 0.05
cfg_stats.clusteralpha     = 0.05; %since we don't correct the tail, later only clusters <0.025 are sig
cfg_stats.numrandomization = 1000;  

Nsub = numel(task_allsubj.timelock.ERP.target_low);
cfg_stats.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg_stats.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg_stats.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg_stats.uvar                = 2; % the 2nd row in cfg.design contains the subject number
cfg_stats.minnbchan  = 2;


%% HEP statistics
cfg_stats.latency     = [-0.2 0.6];
cfg_stats.avgovertime = 'no'; % not average in the ROI

% High vs. low target HEP
stat_HEP_target = ft_timelockstatistics(cfg_stats, task_allsubj.timelock.HEP.target_low{:}, task_allsubj.timelock.HEP.target_high{:});

% High vs. low target HEP correct random pseudotrial
stat_HEP_target_correct_random = ft_timelockstatistics(cfg_stats, task_allsubj.timelock.HEP.target_low_correct_random{:}, task_allsubj.timelock.HEP.target_high_correct_random{:});

% High vs. low target HEP random pseudo
stat_HEP_target_random_pseudo = ft_timelockstatistics(cfg_stats, task_allsubj.timelock.HEP.target_low_random_pseudo{:}, task_allsubj.timelock.HEP.target_high_random_pseudo{:});

% High vs. low target RT HEP random pseudo
stat_HEP_target_RT_random_pseudo = ft_timelockstatistics(cfg_stats, task_allsubj.timelock.HEP.target_low_RT_random_pseudo{:}, task_allsubj.timelock.HEP.target_high_RT_random_pseudo{:});

% High vs. low target HEP RT
stat_HEP_target_RT = ft_timelockstatistics(cfg_stats, task_allsubj.timelock.HEP.target_low_RT{:}, task_allsubj.timelock.HEP.target_high_RT{:});

% High vs. low target HEP RT correct random
stat_HEP_target_RT_correct_random = ft_timelockstatistics(cfg_stats, task_allsubj.timelock.HEP.target_low_RT_correct_random{:}, task_allsubj.timelock.HEP.target_high_RT_correct_random{:});

%% add masks to data
% HEP
task_allsubj.grandAvg.HEP.target_low.mask = stat_HEP_target.mask;
task_allsubj.grandAvg.HEP.target_high.mask = stat_HEP_target.mask;
task_allsubj.grandAvg.HEP.target_low_RT.mask = stat_HEP_target_RT.mask;
task_allsubj.grandAvg.HEP.target_high_RT.mask = stat_HEP_target_RT.mask;
task_allsubj.grandAvg.HEP.target_high_RT_correct_random.mask = stat_HEP_target_RT_correct_random.mask;
task_allsubj.grandAvg.HEP.target_high_RT_correct_random.mask = stat_HEP_target_RT_correct_random.mask;
task_allsubj.grandAvg.HEP.target_high_correct_random.mask = stat_HEP_target_correct_random.mask;
task_allsubj.grandAvg.HEP.target_high_correct_random.mask = stat_HEP_target_correct_random.mask;
task_allsubj.grandAvg.HEP.target_low_RT_random_pseudo.mask = stat_HEP_target_RT_random_pseudo.mask;
task_allsubj.grandAvg.HEP.target_high_RT_random_pseudo.mask = stat_HEP_target_RT_random_pseudo.mask;
task_allsubj.grandAvg.HEP.target_low_random_pseudo.mask = stat_HEP_target_random_pseudo.mask;
task_allsubj.grandAvg.HEP.target_high_random_pseudo.mask = stat_HEP_target_random_pseudo.mask;

%% Save output
save([savepath, 'stats.mat'], 'stat_HEP_target',  ...
    'stat_HEP_target_correct_random', 'stat_HEP_target_random_pseudo', ...
    'stat_HEP_target_RT_random_pseudo', 'stat_HEP_target_RT', 'stat_HEP_target_RT_correct_random');

task_allsubj_grandAvg.grandAvg = task_allsubj.grandAvg;
task_allsubj_grandAvg.ID = task_allsubj.ID;
save([savepath, 'task_grandAvg.mat'], 'task_allsubj_grandAvg');

%% create contrasts
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';

diff_HEP_target = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low, task_allsubj.grandAvg.HEP.target_high);
diff_HEP_random_pseudo = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low_random_pseudo, task_allsubj.grandAvg.HEP.target_high_random_pseudo);
diff_HEP_target_correct_random = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low_correct_random, task_allsubj.grandAvg.HEP.target_high_correct_random);
diff_HEP_target_RT_correct_random = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low_RT_correct_random, task_allsubj.grandAvg.HEP.target_high_RT_correct_random);
diff_HEP_target_RT = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low_RT, task_allsubj.grandAvg.HEP.target_high_RT);
diff_HEP_target_random_pseudo = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low_random_pseudo, task_allsubj.grandAvg.HEP.target_high_random_pseudo);
diff_HEP_target_random_pseudo_RT = ft_math(cfg, task_allsubj.grandAvg.HEP.target_low_RT_random_pseudo, task_allsubj.grandAvg.HEP.target_high_RT_correct_random);

%% Get clusters and their respective averages

selchan = stat_HEP_target_correct_random.label;

cfg_topo = [];
cfg_topo.channel = 1:31; %{'all','-ECG'}; % exclude ECG channel
cfg_topo.highlight = 'on';
cfg_topo.comment   = 'no';
cfg_topo.selchan = stat_HEP_target_correct_random.label;
cfg_topo.layout=layout;
cfg_topo.interactive = 'no';

% List of all comparisons
comparisons = {
    'diff_HEP_random_pseudo', 'task_allsubj.grandAvg.HEP.target_low_random_pseudo', 'task_allsubj.grandAvg.HEP.target_high_random_pseudo', 'stat_HEP_target_random_pseudo';
    'diff_HEP_target_correct_random', 'task_allsubj.grandAvg.HEP.target_low_correct_random', 'task_allsubj.grandAvg.HEP.target_high_correct_random', 'stat_HEP_target_correct_random';
    'diff_HEP_target_RT', 'task_allsubj.grandAvg.HEP.target_low_RT', 'task_allsubj.grandAvg.HEP.target_high_RT', 'stat_HEP_target_RT';
    'diff_HEP_target', 'task_allsubj.grandAvg.HEP.target_low', 'task_allsubj.grandAvg.HEP.target_high', 'stat_HEP_target';
    'diff_HEP_target_RT_correct_random', 'task_allsubj.grandAvg.HEP.target_low_RT_correct_random', 'task_allsubj.grandAvg.HEP.target_high_RT_correct_random', 'stat_HEP_target_RT_correct_random';
    'diff_HEP_target_random_pseudo', 'task_allsubj.grandAvg.HEP.target_low_random_pseudo', 'task_allsubj.grandAvg.HEP.target_high_random_pseudo', 'stat_HEP_target_random_pseudo';
    'diff_HEP_target_random_pseudo_RT', 'task_allsubj.grandAvg.HEP.target_low_RT_random_pseudo', 'task_allsubj.grandAvg.HEP.target_high_RT_random_pseudo', 'stat_HEP_target_RT_random_pseudo'
    };

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';

for i = 1:size(comparisons, 1)
  
    % Compute the difference
    diff_name = comparisons{i, 1};
    cond1 = eval(comparisons{i, 2});
    cond2 = eval(comparisons{i, 3});
    stat_name = comparisons{i, 4};

    diff_data = ft_math(cfg, cond1, cond2);
    eval([diff_name ' = diff_data;']);

    if ~isfield(eval(stat_name), 'posclusters') && isfield(eval(stat_name), 'negclusters')
        continue
    end

    % Plot significant topographies
    direction = {'pos', 'neg'};
    for dir = direction
        direct = dir{1};

        probSum = 0;

        switch direct
            case 'pos'
                if ~isfield(eval(stat_name), 'posclusters')
                    continue
                end
                fullVarName = [stat_name '.posclusters.prob'];
                probSum = sum(evalin('base', ['[', fullVarName, ']'])<= 0.025);
            case 'neg'
                if ~isfield(eval(stat_name), 'negclusters')
                    continue
                end
                fullVarName = [stat_name '.negclusters.prob'];
                probSum = sum(evalin('base', ['[', fullVarName, ']'])<= 0.025);
        end

        if probSum == 0
            continue
        end

        for sig = 1:probSum
            cfg_topo = [];
            cfg_topo.channel = 1:31; % Exclude ECG channel
            cfg_topo.highlight = 'on';
            cfg_topo.comment   = 'no';
            cfg_topo.selchan = eval([stat_name '.label']);
            cfg_topo.layout = layout;
            cfg_topo.interactive = 'no';


            cfg_topo = plot_sig_topo(eval(stat_name), cfg_topo, diff_data, diff_data.time, sig, direct, 1);
            print('-dsvg', [savepath, diff_name, '_topo_', num2str(sig), '_', direct, '.svg']);

            plot_erp_var(cond1, cond2, [-1.52 1.5], cfg_topo.highlightchannel, [], 'sem', eval(stat_name), direct, ...
                {[diff_name, ' task'], [diff_name, ' rest']}, 'northwest', [], savepath, [diff_name, '_', num2str(sig), '_', direct, '.svg']);
        end
    end

    %save cluster information as table
    save_cluster_info(eval(stat_name), stat_name, savepath)
    close all
end


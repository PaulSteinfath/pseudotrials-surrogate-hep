%% MAIN_PREPROCESSING_PIPELINE - Oddball HEP Analysis
%
% This script implements the main preprocessing pipeline for analyzing EEG
% data in an Oddball HEP (Heartbeat Evoked Potential) experiment. The
% pipeline includes:
%   - Path initialization and directory setup.
%   - Preprocessing of both task and resting-state EEG data.
%   - Initial preprocessing steps: filtering, downsampling, and ICA.
%   - ICA component selection and application to filtered data.
%   - Creating pseudo trials for task and rest conditions.
%   - Statistical analysis on preprocessed data.
%
% Key steps:
% 1. Path setup and library initialization (e.g., EEGLAB, HEPLAB, Fieldtrip).
% 2. File management: ensuring the required directories exist.
% 3. Preprocessing for task and rest EEG data.
% 4. Artifact correction using ICA.
% 5. Sorting trials and statistical analysis.
%
% Requirements:
%   - EEGLAB toolbox with bva_io, erp_lab plugins, and HEPLAB and Fieldtrip dependencies.
%   - EEG data in BrainVision format (.vhdr, .vmrk, .eeg).
%
% Outputs:
%   - Preprocessed EEG data (filtered, ICA-applied) saved in project directories.
%   - Logfiles for errors encountered during preprocessing.
%   - Pseudo trials and statistical results.
%
% Author:
%   Paul Steinfath
%
% Example:
%   Run the script to process task and rest EEG data:
%     1. Filter and preprocess data.
%     2. Perform ICA and correct artifacts.
%     3. Sort trials and perform statistical analysis.

%% Initialize workspace
clc; clear all; close all; % Clear command window, variables, and figures.

%% Initialize paths
% Base paths
eeglab_path = '';
code_path = '';


% Add required paths
addpath(eeglab_path);
addpath([code_path, 'functions']);
addpath([code_path, 'Preprocessing/']);
addpath([code_path, 'Stats/Oddball/']);
addpath([code_path, 'Simulations/']);

addpath('/HEPLAB-master/HEPLAB-master/Functions'); % HEPLAB path
addpath('/fieldtrip-20230822/'); % Fieldtrip path
addpath(genpath([code_path, 'functions/boundedline'])); % Bounded lines path

% Raw data paths
task_raw_path = ''; %oddball task data


rest_raw_path = {'',...
    '',...
    ''};

% Project paths
% Task
task_basedir = '/final/Task/';
task_pre_ica_path = [task_basedir, 'ICA/'];
task_post_ica_path = [task_basedir, 'postICA/'];
task_filtered_path = [task_basedir, 'postICA_03_45Hz/'];
task_epochs_path = [task_basedir, 'epochs/'];
task_output_path = [task_basedir, 'output/'];
task_qa_path = [task_basedir, 'QA/'];
task_error_path = [task_basedir, 'Logfiles/'];
settings_path = [code_path, 'settings/']; % same settings for rest / task data

% Rest
rest_basedir = '/final/Rest/';
rest_pre_ica_path = [rest_basedir, 'ICA/'];
rest_post_ica_path = [rest_basedir, 'postICA/'];
rest_filtered_path = [rest_basedir, 'postICA_03_45Hz_markers_added/'];
rest_epochs_path = [rest_basedir, 'epochs/'];
rest_output_path = [rest_basedir, 'output/'];
rest_qa_path = [rest_basedir, 'QA/'];
rest_error_path = [rest_basedir, 'Logfiles/'];

% Task vs Rest
task_vs_rest_path = '/final/task_vs_rest';

% Simulation
sim_basedir = '/final/simulations/permutations/';

% Create all task directories
dirs_to_create = {task_pre_ica_path, task_post_ica_path, task_filtered_path, ...
    task_epochs_path, task_qa_path, task_error_path};
create_dirs(dirs_to_create);

% Create rest directories
rest_dirs_to_create = {rest_pre_ica_path, rest_post_ica_path, rest_filtered_path, ...
    rest_qa_path, rest_epochs_path, rest_error_path, task_vs_rest_path};
create_dirs(rest_dirs_to_create);

% Create simulation directories
create_dirs({sim_basedir});

% Initialize EEGLAB
eeglab; close;

%% Preprocessing settings
% Down-sampling frequency in Hz
fs = 250;

% Location of standard electrode positions
elecfile = [eeglab_path, '/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp'];
% High- and Low-pass filtering cutoffs, in Hz
highpass_cu = 0.3;
lowpass_cu = 45;
% Criterion for cleaning channel, flat and line-noise 
chan_crit = 0.85; % correlation with neighboring electrodes
flatline_crit = 5; % in s
ln_crit = 3;

%% ICA cleaning settings
% Filtering for ICA decomposition, in Hz
ica_highpass_cu = 1;
ica_lowpass_cu = lowpass_cu;
% Notch filtering
line_noise_f = 50; % line noise, in Hz
% ICA and ECG epoching for ECG-related artifact detection, in s
ica_window = [-0.050 0.600];

% SD threshold for ECG-related artifact detection
sd_ecg_thresh = 1.5; % ica_thresh
% threshold for correlation between ECG and ICA components
ecg_tresh = 0.8;
% threshold probability for IClabeling - muscle
muscle_thresh = 0.5;
% threshold probability for IClabeling - eye
eye_thresh = 0.6;
% threshold probability for IClabeling - line noise
ln_thresh = 0.5;
% threshold probability for IClabeling - channel noise
chann_thresh = 0.4;
% threshold probability for IClabeling - other
other_thresh = 0.5;

% Create a dictionary using containers.Map
thresholds = containers.Map(...
    {'sd_ecg', 'ecg', 'muscle', 'eye', ...
     'ln', 'chann', 'other'}, ...
    {sd_ecg_thresh, ecg_tresh, muscle_thresh, eye_thresh, ...
     ln_thresh, chann_thresh, other_thresh});

%% Pseudo-trials settings
% Binlister including new random pseudotrials
binlister_file = [settings_path, 'Oddball-LIFE-binlister_BK.txt'];
% Epoch length
epoch_length = [-1500, 1500];
% LDA matrices
lda_file = [settings_path, 'LDA/erp_t_filter.mat'];
ids_file = [settings_path, 'LDA/ids.txt'];
% inter-stimulus-interval, in s
isi = 1.5; 
% minimum inter-stimulus-interval that is allowed, events that are faster will be removed
min_isi = 0.5; % in s
% select only HEP occuring in hep_time_window BEFORE stimulus
hep_prestim_window = [1.1, 0.6]; % in s
% HEP epoching
hep_epoch = [250 600]; % in ms
% Threshold noisy epoch
noise_epoch_tresh = 150; % in uV
% Baseline correction
rm_base = 1; % yes = 1, no = 0
baseline = [-150 -50]; % in ms

%% Markers settings
% task, response and target stimulus markers in recording
target_marker = '20';
before_target_marker = '22';
response_marker = {'S  1'};
standard_stim_marker = {'S 10'};
target_stim_marker = {['S ' target_marker]};
novelty_stim_marker = {'S 30'};
% R marker in recording
r_marker = '88';
% pseudo trial marker
pseudo_marker = '4';
% random trial marker
random_marker = '8';
% marker for trials with low amplitude target
low_target_marker = '1';
% marker for trials with high amplitude target
high_target_marker = '2';
% marker for low reaction time
low_rt_marker = '71';
% marker for high reaction time
high_rt_marker = '72';
% marker for HEP_target_low_RT_random_pseudo
HEP_target_low_RT_random_pseudo_marker = '15';
% marker for HEP_target_low_RT_random_pseudo
HEP_target_high_RT_random_pseudo_marker =  '16';
% marker for HEP_target_low_RT_random_pseudo
ERP_target_low_RT_pseudo_marker =  '17';
% marker for HEP_target_low_RT_random_pseudo
ERP_target_high_RT_pseudo_marker =  '18';
% markers for permutation shuffling
permutations_markers = {'3', '4'};

% Dictionary of markers
markers = containers.Map(...
    {'target', 'before_target', 'response', ...
     'standard_stim', 'target_stim', 'novelty_stim', 'r', 'pseudo', ...
     'random', 'low_target', 'high_target', ...
     'low_rt', 'high_rt', ...
     'HEP_target_low_RT_random_pseudo', ...
     'HEP_target_high_RT_random_pseudo', ...
     'ERP_target_low_RT_random_pseudo', ...
     'ERP_target_high_RT_random_pseudo', ...
     'permutations'}, ...
    {target_marker, before_target_marker, response_marker, ...
     standard_stim_marker,target_stim_marker, novelty_stim_marker, r_marker, pseudo_marker, ...
     random_marker, low_target_marker, high_target_marker, ...
     low_rt_marker, high_rt_marker, ...
     HEP_target_low_RT_random_pseudo_marker, ...
     HEP_target_high_RT_random_pseudo_marker, ...
     ERP_target_low_RT_pseudo_marker, ...
     ERP_target_high_RT_pseudo_marker, ...
     permutations_markers});

%% Simulation settings

% Binlister file
sim_binlister_file = [settings_path 'sim_perm_BINLISTER.txt'];

% Interval between stimuli, in seconds
stim_interval = 4;

% Simulation parameters
sim_n_subjects = 50;
sim_n_permutations = 50;
sim_relationship_types = {'direct', 'inverse'};
sim_relationship_scales = [1, 0.8, 0.6, 0.4, 0.2, 0];

% Resampling frequency, in Hz
sim_fs = 250;

% Baseline correction
sim_rm_base = 1;
sim_baseline = [-150 -50]; % in ms

% Statistical settings
sim_alpha = 0.05;
sim_cluster_alpha = 0.05;
sim_num_randomization = 1000;

% Channel and timing settings
sim_channel = 'Pz';
sim_latency = [-0.2 0.6];
sim_epoch_length = [-1500, 1500];
sim_hep_time_window = [1.1, 0.6];

% Parameters for simulating responses with phase-randomized ICA data
sim_params = struct();
sim_params.triggers      = markers;
sim_params.hep_params    = [1.5, 0.8];
sim_params.p300_params   = [7, 2];
sim_params.p300_timing   = [0.2, 0.6];
sim_params.hep_timing    = [0.2, 0.4];


%% Step 1: Initial preprocessing and ICA - Task
fprintf('Running step 1: Initial preprocessing and ICA for Task\n');
a_1_oddball_preprocessing(task_raw_path, task_pre_ica_path, task_error_path, fs, elecfile, ica_highpass_cu, ica_lowpass_cu, line_noise_f, flatline_crit)

%% Step 1b: Initial preprocessing and ICA - Rest
fprintf('Running step 1b: Initial preprocessing and ICA for Rest\n');
a_1_oddball_preprocessing(rest_raw_path, rest_pre_ica_path, rest_error_path, fs, elecfile, ica_highpass_cu, ica_lowpass_cu, line_noise_f, flatline_crit)

%% Step 2: Select ICA components - Task
fprintf('Running step 2: ICA component selection\n');
a_2_select_ICA_components(task_pre_ica_path, task_post_ica_path, task_error_path, task_qa_path, ica_window, thresholds, markers)

%% Step 2b: Select ICA components - Rest
fprintf('Running step 2: ICA component selection for Rest\n');
a_2_select_ICA_components(rest_pre_ica_path, rest_post_ica_path, rest_error_path, rest_qa_path, ica_window, thresholds, markers)

%% Step 3: Apply ICA to filtered data - Task
fprintf('Running step 3: Applying ICA to filtered data\n');
a_3_apply_ICA_components_to_filtered_data(task_raw_path, task_filtered_path, task_pre_ica_path, task_post_ica_path, task_error_path, task_qa_path, fs, elecfile, highpass_cu, lowpass_cu, line_noise_f, chan_crit, ln_crit)


%% Step 3b: Apply ICA to filtered data - Rest
fprintf('Running step 3: Applying ICA to filtered rest data\n');
a_3b_apply_ICA_components_to_filtered_data_add_task_marker(task_raw_path, rest_raw_path, rest_filtered_path, rest_pre_ica_path, rest_post_ica_path, rest_error_path, rest_qa_path, fs, elecfile, highpass_cu, lowpass_cu, line_noise_f)


%% Step 4: Sort and create pseudo trials - Task
fprintf('Running step 4: Sorting and creating pseudo trials\n');
a_4_sort_random_pseudo(task_filtered_path, task_epochs_path, task_error_path, task_qa_path, binlister_file, epoch_length,lda_file, ids_file, markers, min_isi, isi, hep_prestim_window, hep_epoch, noise_epoch_tresh)


%% Step 4b: Sort and create pseudo trials - Rest
fprintf('Running step 4: Sorting and creating pseudo trials for Rest\n');
a_4_sort_random_pseudo(rest_filtered_path, rest_epochs_path, rest_error_path, rest_qa_path, binlister_file, epoch_length,lda_file, ids_file, markers, min_isi, isi, hep_prestim_window, hep_epoch, noise_epoch_tresh)

%% Step 5: Run stats - Task
fprintf('Running step 5: Stats\n');
task = 1; % 0: rest data, 1: task data
run([code_path, 'Stats/Oddball/stats_sort_random_pseudo.m']);


%% Step 5: Run stats - Rest
fprintf('Running step 5: Stats for Rest\n');
task = 0; % 0: rest data, 1: task data
run([code_path, 'Stats/Oddball/stats_sort_random_pseudo.m']);


%% Step 6: Run simulation and permutation analysis
fprintf('Running simulation and permutation analysis\n');
run_permutations_phase_rand(task_filtered_path, sim_basedir, sim_n_subjects, sim_n_permutations, ...
    sim_relationship_types, sim_relationship_scales, sim_rm_base, sim_baseline, sim_fs, sim_alpha, ...
    sim_cluster_alpha, sim_num_randomization, sim_channel, sim_latency, ...
    sim_epoch_length, stim_interval, sim_hep_time_window, sim_params, settings_path, hep_epoch, lda_file, ids_file, binlister_file, markers);

%% Step 7: Run permutation analysis on preprocessed files
fprintf('Running permutation analysis on preprocessed files\n');

% Number of permutations to run - 10 for testing
% num_perms = 1000; % for real analysis
num_perms = 10; 

% Define base directory for permutation analysis rest & task
perm_basedir_task = fullfile(task_basedir, 'permutations');
perm_basedir_rest = fullfile(rest_basedir, 'permutations');
create_dirs({perm_basedir_task, perm_basedir_rest});

% Path to the binlister file for permutation analysis - smaller than main binlister file for processing speed
perm_binlister = [settings_path, 'Oddball-LIFE-binlister_perm.txt'];

% Run the permutation analysis - Task
surrogate_statistics(task_filtered_path, perm_basedir_task, task_error_path,  ...
    perm_binlister, epoch_length, lda_file, ids_file, markers, min_isi, isi, ...
    hep_prestim_window, hep_epoch, noise_epoch_tresh, rm_base, baseline, num_perms, settings_path);
       
% Run the permutation analysis - Rest
surrogate_statistics(rest_filtered_path, perm_basedir_rest, rest_error_path,  ...
    perm_binlister, epoch_length, lda_file, ids_file, markers, min_isi, isi, ...
    hep_prestim_window, hep_epoch, noise_epoch_tresh, rm_base, baseline, num_perms, settings_path);

%% Step 8: Summarize Permutation Statistics
% Define directories for stat files and image save path
perm_img_save_path = fullfile(perm_basedir_task, 'stats_img');

% Call the summarize_stats_abs function
summarize_stats_abs(perm_basedir_task, perm_img_save_path, task_output_path, rm_base);

fprintf('Permutation statistics summarized!\n');

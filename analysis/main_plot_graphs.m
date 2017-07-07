clear; close all; 

%% data directory to plot
id = '06141205';
% id = '06150910';

%% some parameters

% total number of files; number of time points in each file; duration per file
% K = 2000; T = 50; T0 = 0.1;
K = 200; T = 500; T0 = 1;

% for coherence split the entire course into P periods
P = 20; per_P = 5000;
% P = 39; per_P = 5000; % each period is 10s, overlapping 5s with the previous period

%%  input directory and file names

folder_list = dir('./data');
for i = 1 : length(folder_list)
    if ~isempty(strfind(folder_list(i).name, id))
        DATA_DIR = ['./data/' folder_list(i).name '/'];
        
        if ~isempty(strfind(folder_list(i).name, 'brain'))
            type = 'brain';
        else
            type = 'sphere';
        end
    end
end
if ~exist('DATA_DIR', 'var')
    error('folder ID does not exist');
end

if exist([DATA_DIR 'raw/'], 'dir')
    RAW_DIR = [DATA_DIR 'raw/'];
else
    RAW_DIR = DATA_DIR;
end
META_FILE = [DATA_DIR 'meta.mat'];

%% load vertices that are closest to electrodes
if strcmp(type, 'sphere')
    [focus_idx, macro_pos, macro_transform, macro_2d, ...
                micro_pos, micro_transform, micro_2d] = ...
        generate_electrode_grid_sphere(RAW_DIR);
else
    loc_grid_center = [64.41, -7.28, 21.48];        % center of the ECoG grid (mm)
    dist_grid = 12;                                 % distance between electrodes (mm)
    find_electrode_grid_center(RAW_DIR, 12);
    flag_dipole = true;
    keyboard;
    
    [focus_idx, macro_pos, macro_transform, macro_2d, ...
                micro_pos, micro_transform, micro_2d] = ...
        generate_electrode_grid_brain(loc_grid_center, dist_grid, RAW_DIR, flag_dipole);
end

%% output directory and file names
if strcmp(type, 'sphere')
    ANALYSIS_DIR = DATA_DIR;
else
    NOTE = 'test_close';
    ANALYSIS_DIR = [DATA_DIR 'grid_'...
                             num2str(loc_grid_center(1), '%.2f') '_' ...
                             num2str(loc_grid_center(2), '%.2f') '_' ...
                             num2str(loc_grid_center(3), '%.2f') '_' num2str(dist_grid) '_' NOTE '/'];
    mkdir(ANALYSIS_DIR);
end

ELEC_FILE = [ANALYSIS_DIR 'electrode_data.mat'];
save(ELEC_FILE, 'focus_idx', 'macro_pos', 'macro_transform', 'macro_2d', ...
                             'micro_pos', 'micro_transform', 'micro_2d');

SAMPLE_DATA_FILE = [ANALYSIS_DIR 'sample_data.mat'];
COHERENCE_FILE = [ANALYSIS_DIR 'coherence.mat'];
VIDEO_FILE = [ANALYSIS_DIR 'movie_sparse_sphere.mp4'];
TRACES_FIG = [ANALYSIS_DIR 'traces.fig'];
COHERENCE_FIG = [ANALYSIS_DIR 'coherence_summary.fig'];

%% load or generate data
if exist(SAMPLE_DATA_FILE, 'file') == 2
    load(SAMPLE_DATA_FILE);
    load(ELEC_FILE);
else
    prepare_sample_data;
end

%% load or compute coherence data
if exist(COHERENCE_FILE, 'file') == 2
    load(COHERENCE_FILE);
else
    prepare_coherence;
end

%% plot seizure course
if strcmp(type, 'sphere')
    plot_course;
end
    
%% plot firing rate and voltage traces
plot_traces;

%% plot coherence statistics
[t_coh,t_coh_conf,t_phi] = deal(macro_t_coh, macro_t_coh_conf, macro_t_phi);
% [t_coh,t_coh_conf,t_phi] = deal(micro_t_coh, micro_t_coh_conf, micro_t_phi);
central_t = 5:10:195;
period_idx = find(central_t == 115); %if drawing one period, draw this
plot_coherence;
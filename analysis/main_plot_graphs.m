clear; close all; 

%% data directory to plot
% id = '07191052';
id = '07181403';
% id = '06141205';
% id = '06150910';

%% some parameters

% total number of files; number of time points in each file; duration per file
K = 2000; T = 50; T0 = 0.1;
% K = 150; T = 500; T0 = 1;

total_time = K*T0;
sparse_time = (1:K) * T0;
fine_time = (1:K*T) * (T0/T);

% for coherence split the entire course into P periods
% P = 100; per_P = 1000;
P = 20; per_P = 5000;
% P = 39; per_P = 5000; % each period is 10s, overlapping 5s with the previous period

%% input directory and file names

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
    NOTE = 'closest7_avg';
    loc_grid_center = [64.41, -7.28, 21.48];        % center of the ECoG grid (mm)
    dist_grid = 12;                                 % distance between electrodes (mm)
    find_electrode_grid_center(RAW_DIR, dist_grid);
    flag_dipole = false;
    closest_N = 7;
    keyboard;
    
    [focus_idx, macro_pos, macro_transform, macro_2d, ...
                micro_pos, micro_transform, micro_2d] = ...
        generate_electrode_grid_brain(loc_grid_center, dist_grid, RAW_DIR, ...
                                      flag_dipole, closest_N);
end

%% output directory and file names
if strcmp(type, 'sphere')
    ANALYSIS_DIR = DATA_DIR;
else
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
VIDEO_FILE = [ANALYSIS_DIR 'movie.mp4'];

COURSE_FIG = [ANALYSIS_DIR 'course.fig'];
TRACES_FIG = [ANALYSIS_DIR 'traces.fig'];
SINGLE_FIG = [ANALYSIS_DIR 'single.fig'];
COHERENCE_FIG = [ANALYSIS_DIR 'coherence_summary.fig'];

%% load grid
if strcmp(type, 'sphere')
    load('N10242_R10.mat');
    pos_hemi = locs(:,3) >= 0;
    neg_hemi = locs(:,3) < 0;
else
    load('N40962.mat');
    load('unitsphere.mat');
    
    surf.vertices = locs;
    surf.faces = tri;
    surf_sphere.vertices = 10 * coord';
    surf_sphere.faces = tri;
end

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
plot_course;

%% plot firing rate and voltage traces
plot_traces;

%% plot single node dynamics
fg = figure;
plot(sparse_time, table2array(single_node));
legend('Qe', 'Qi', 'Ve', 'Vi', 'D22', 'dVe', 'dVi', 'K');
saveas(fg, SINGLE_FIG);

%% plot coherence statistics
[t_coh,t_coh_conf,t_phi] = deal(macro_t_coh, macro_t_coh_conf, macro_t_phi);
% [t_coh,t_coh_conf,t_phi] = deal(micro_t_coh, micro_t_coh_conf, micro_t_phi);
central_t = int32(total_time * (1/P/2 : 1/P : 1-1/P/2));
period_idx = find(central_t == 115); %if drawing one period, draw this
plot_coherence;
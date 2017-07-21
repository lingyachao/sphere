clear;
% clearvars -except id K T0;
close all;
% set(0,'DefaultFigureWindowStyle','docked');

%% data directory to plot
if ~exist('id', 'var')
    id = '07191525';
    % id = '07191052';
    % id = '07181403';
    % id = '06141205';
    % id = '06150910';
end

% brain type only
NOTE = 'closest7_avg';
loc_grid_center = [64.41, -7.28, 21.48];        % center of the ECoG grid (mm)
dist_grid = 12;                                 % distance between electrodes (mm)

%% *** SPECIFY *** input directory and file names
folder_list = dir('./data');
for i = 1 : length(folder_list)
    folder_name = folder_list(i).name;
    
    if ~isempty(strfind(folder_name, id))
        DATA_DIR = ['./data/' folder_name '/'];
        
        if ~isempty(strfind(folder_name, 'brain'))
            type = 'brain';
        else
            type = 'sphere';
        end
    end
end

RAW_DIR = [DATA_DIR 'raw/'];
META_FILE = [DATA_DIR 'vars.mat'];
load(META_FILE);

%% *** SPECIFY *** output directory and file names
if strcmp(type, 'sphere')
    ANALYSIS_DIR = DATA_DIR;
else
    ANALYSIS_DIR = [DATA_DIR 'grid_'...
                             num2str(loc_grid_center(1), '%.2f') '_' ...
                             num2str(loc_grid_center(2), '%.2f') '_' ...
                             num2str(loc_grid_center(3), '%.2f') '_' ...
                             num2str(dist_grid) '_' NOTE '/'];
    mkdir(ANALYSIS_DIR);
end

VIDEO_FILE = [ANALYSIS_DIR 'movie.mp4'];

SAMPLE_DATA_FILE = [ANALYSIS_DIR 'data_sample.mat'];
COHERENCE_FILE   = [ANALYSIS_DIR 'data_coherence.mat'];
ELEC_FILE        = [ANALYSIS_DIR 'data_electrode.mat'];

COURSE_FIG       = [ANALYSIS_DIR 'fig_course.fig'];
TRACES_FIG       = [ANALYSIS_DIR 'fig_traces.fig'];
SINGLE_FIG       = [ANALYSIS_DIR 'fig_single.fig'];
COH_MACRO_FIG    = [ANALYSIS_DIR 'fig_coh_summary_macro.fig'];
COH_MICRO_FIG    = [ANALYSIS_DIR 'fig_coh_summary_micro.fig'];

%% *** SPECIFY *** time sequences
K = length(dir([RAW_DIR 'seizing_*.mat']));
T0 = 1 * (K<1000) + 0.1 * (K >= 1000);
T = 500 * T0;

total_time = K*T0;
sparse_time = (1:K) * T0;
fine_time = (1:K*T) * (T0/T);

% for coherence split the entire course into P periods
P = total_time/10;
per_P = K*T/P;

%% *** SPECIFY *** grid
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

%% *** LOAD/GENERATE *** sample data
if exist(SAMPLE_DATA_FILE, 'file') == 2
    load(ELEC_FILE);
    load(SAMPLE_DATA_FILE);
else
    prepare_electrode_data;
    prepare_sample_data;
end

%% *** LOAD/GENERATE *** coherence data
if exist(COHERENCE_FILE, 'file') == 2
    load(COHERENCE_FILE);
else
    prepare_coherence;
end

%% *** PLOT *** seizure course
plot_course;

%% *** PLOT *** firing rate and voltage traces
plot_traces;

%% *** PLOT *** single node dynamics
fg_single = figure;
plot(sparse_time, table2array(single_node));
legend('Qe', 'Qi', 'Ve', 'Vi', 'D22', 'dVe', 'dVi', 'K');
saveas(fg_single, SINGLE_FIG);

%% *** PLOT *** coherence statistics
central_t = int32(total_time * (1/P/2 : 1/P : 1-1/P/2));
period_idx = length(central_t) - 2; % estimate wave for this period

COHERENCE_FIG = COH_MACRO_FIG;
[t_coh,t_coh_conf,t_phi,electrode_2d] = deal( ...
    macro_t_coh, macro_t_coh_conf, macro_t_phi, macro_2d);
fprintf('macro electrodes ');
plot_coherence;

if ~isempty(micro_pos)
    COHERENCE_FIG = COH_MICRO_FIG;
    [t_coh,t_coh_conf,t_phi,electrode_2d] = deal( ...
        micro_t_coh, micro_t_coh_conf, micro_t_phi, micro_2d);
    fprintf('micro electrodes ');
    plot_coherence;
end

%% *** CALCULATE *** recruitment speed
[~, node_e] = min(macro_pos(:,3));
[~, node_l] = max(macro_pos(:,3));

t_e = fine_time(find(Qe_macro(:,node_e) > 15, 1));
t_l = fine_time(find(Qe_macro(:,node_l) > 15, 1));

travel_dist = dist(macro_2d([node_e, node_l], :)');
recruitment_speed = travel_dist(1,2) / (t_l-t_e);
fprintf(['recruitment speed is ' num2str(recruitment_speed) ' cm/s \n']);

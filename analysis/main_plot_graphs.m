clearvars -except id K T0;
close all; 

%% data directory to plot
if ~exist('id', 'var')
    id = '07191525';
    % id = '07191052';
    % id = '07181403';
    % id = '06141205';
    % id = '06150910';
end

%% input directory and file names
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

%% time sequences
K = length(dir([RAW_DIR 'seizing_*.mat']));
T0 = 1 * (K<1000) + 0.1 * (K >= 1000);
T = 500 * T0;

total_time = K*T0;
sparse_time = (1:K) * T0;
fine_time = (1:K*T) * (T0/T);

% for coherence split the entire course into P periods
P = 20; per_P = K*T/P;

%% gain matrix for electrodes
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

VIDEO_FILE = [ANALYSIS_DIR 'movie.mp4'];

SAMPLE_DATA_FILE = [ANALYSIS_DIR 'data_sample.mat'];
COHERENCE_FILE   = [ANALYSIS_DIR 'data_coherence.mat'];
ELEC_FILE        = [ANALYSIS_DIR 'data_electrode.mat'];
save(ELEC_FILE, 'focus_idx', 'macro_pos', 'macro_transform', 'macro_2d', ...
                             'micro_pos', 'micro_transform', 'micro_2d');

COURSE_FIG       = [ANALYSIS_DIR 'fig_course.fig'];
TRACES_FIG       = [ANALYSIS_DIR 'fig_traces.fig'];
SINGLE_FIG       = [ANALYSIS_DIR 'fig_single.fig'];
COH_MACRO_FIG    = [ANALYSIS_DIR 'fig_coh_summary_macro.fig'];
COH_MICRO_FIG    = [ANALYSIS_DIR 'fig_coh_summary_micro.fig'];

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
central_t = int32(total_time * (1/P/2 : 1/P : 1-1/P/2));
period_idx = find(central_t == 175); %if drawing one period, draw this

COHERENCE_FIG = COH_MACRO_FIG;
[t_coh,t_coh_conf,t_phi,electrode_2d] = deal( ...
    macro_t_coh, macro_t_coh_conf, macro_t_phi, macro_2d);
plot_coherence;

if ~isempty(micro_pos)
    COHERENCE_FIG = COH_MICRO_FIG;
    [t_coh,t_coh_conf,t_phi,electrode_2d] = deal( ...
        micro_t_coh, micro_t_coh_conf, micro_t_phi, micro_2d);
    plot_coherence;
end

%% calculate recruitment speed
[~, node_e] = min(macro_pos(:,3));
[~, node_l] = max(macro_pos(:,3));

t_e = fine_time(find(Qe_macro(:,node_e) > 15, 1));
t_l = fine_time(find(Qe_macro(:,node_l) > 15, 1));

travel_dist = dist(macro_2d([node_e, node_l], :)');
recruitment_speed = travel_dist(1,2) / (t_l-t_e);
fprintf(['recruitment speed is ' num2str(recruitment_speed) ' cm/s \n']);

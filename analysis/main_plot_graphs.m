clear; close all; 

% data directory to plot
id = '06141205';

% get dir full-name
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

% define subfolder/subfile names
if exist([DATA_DIR 'raw/'], 'dir')
    RAW_DIR = [DATA_DIR 'raw/'];
else
    RAW_DIR = DATA_DIR;
end
META_FILE = [DATA_DIR 'meta.mat'];
ELEC_FILE = [DATA_DIR 'electrode_data.mat'];
SAMPLE_DATA_FILE = [DATA_DIR 'sample_data.mat'];
COHERENCE_FILE = [DATA_DIR 'coherence.mat'];
VIDEO_FILE = [DATA_DIR 'movie_sparse_sphere.mp4'];
TRACES_FIG = [DATA_DIR 'traces.fig'];
COHERENCE_FIG = [DATA_DIR 'coherence_summary.fig'];

% total number of files; number of time points in each file; duration per file
% K = 2000; T = 50; T0 = 0.1;
K = 200; T = 500; T0 = 1;

% for coherence split the entire course into P periods
P = 20; per_P = 5000;
% P = 39; per_P = 5000; % each period is 10s, overlapping 5s with the previous period

%% load or generate data
if exist(SAMPLE_DATA_FILE, 'file') == 2
    load(SAMPLE_DATA_FILE);
    load(ELEC_FILE);
elseif strcmp(type, 'sphere')
    prepare_sample_data;
else
    prepare_sample_data_brain;
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
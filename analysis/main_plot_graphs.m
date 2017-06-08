clear; close all; 

% data directory to plot
id = '06041833';

folder_list = dir('./data');
for i = 1 : length(folder_list)
    if ~isempty(strfind(folder_list(i).name, id))
        DATA_DIR = ['./data/' folder_list(i).name '/'];
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
SAMPLE_DATA_FILE = [DATA_DIR 'sample_data.mat'];
COHERENCE_FILE = [DATA_DIR 'coherence.mat'];
VIDEO_FILE = [DATA_DIR 'movie_sparse_sphere.mp4'];

% add functions to path
addpath(genpath('./analysis'));

% total number of files; number of time points in each file; duration per file
K = 2000; T = 50; T0 = 0.1;
% K = 200; T = 500; T0 = 1;

% for coherence split the entire course into P periods
% P = 39; % each period is 10s, overlapping 5s with the previous period
% per_P = 5000; % size(Qe_macro,1) / P;
P = 20;
per_P = 5000;

% load or generate data
if exist(SAMPLE_DATA_FILE, 'file') == 2
    load(SAMPLE_DATA_FILE);
else
    prepare_sample_data;
    % error('script should be run in the root folder'); 
end

% load or compute coherence data
if exist(COHERENCE_FILE, 'file') == 2
    load(COHERENCE_FILE);
else
    prepare_coherence;
    % error('script should be run in the root folder');
end

% plot seizure course
plot_course;

% plot firing rate and voltage traces
plot_traces;

% plot coherence statistics
% [t_coh,t_coh_conf,t_phi] = deal(macro_t_coh, macro_t_coh_conf, macro_t_phi);
[t_coh,t_coh_conf,t_phi] = deal(micro_t_coh, micro_t_coh_conf, micro_t_phi);
plot_coherence;
clear; close all;

% specify the data set to analyze
data_set = 'sphere_N10242_R10_03301320_high_inhibit_1.2';
DATA_DIR = ['./data/' data_set '/'];

%% downsample the data
DOWNSAMPLE_FILE = [DATA_DIR 'downsample_data.mat'];
if exist(DOWNSAMPLE_FILE, 'file') == 2
    load(DOWNSAMPLE_FILE);
else
    [data_macro,data_micro,fs] = create_sample(DATA_DIR);
end

%% compute average coherence over time

% [avg_coh_macro, avg_coh_micro] = coherence_with_time(data_macro, data_micro, fs);

% % choose micro or macro to analyze
% data = data_macro;
% time = 1/fs * (0 : 1 : size(data,1) - 1);
% 
% y = [-12, -12, -12,   0,  0,  0,  12, 12, 12];
% x = [-12,   0,  12, -12,  0, 12, -12,  0, 12];
% position = [x;y]';
% 
% figure;
% plot(time, data);
% xlabel('Time (s)');
% ylabel('Voltage (uV)');
% 
% figure;
% plot(position(:,1), position(:,2), 'o');
% xlabel('Electrode position (mm)');
% ylabel('Electrodes position (mm)');


%% compute coherence and phase between 120-130s
data = data_macro(:,:,12);

CHRONUX_PATH = './chronux_2_12'; % required toolbox
addpath(genpath(CHRONUX_PATH));

BAND = [1 13];                  % Select a frequency range to analyze
TW = 20;                        % Time-bandwidth product
ntapers = 2*TW-1;               % Choose the # of tapers.
params.tapers = [TW, ntapers];  % ... time-bandwidth product and tapers.
params.Fs = fs;                 % ... sampling rate
params.pad = -1;                % ... no zero padding.
params.fpass = BAND;            % ... freq range to pass
params.err = [1 0.05];          % ... theoretical error bars, p=0.05.

[coh, phi, freq, coh_conf] = compute_coherence(data, params);
save(['./data/' data_set '/coherence'], 'coh', 'phi', 'freq', 'coh_conf');

figure;
plot(freq, squeeze(coh(1,8,:)));
xlabel('Frequency (Hz)');
ylabel('Coherence');

figure;
subplot(1,2,1)
imagesc(squeeze(coh(:,:, find(freq >= 6, 1))));
c = colorbar;
c.Label.String = 'Coherence at 6 Hz';
xlabel('Electrodes');
ylabel('Electrodes');
subplot(1,2,2)
imagesc(squeeze(phi(:,:, find(freq >= 6, 1))));
c = colorbar;
c.Label.String = 'Phase at 6 Hz';
xlabel('Electrodes');
ylabel('Electrodes');


%% estimate delays between electrodes
[delay, delay_ci_lo, delay_ci_up] = compute_delay(coh, coh_conf, phi, freq);

figure;
imagesc(1000 * delay);
c = colorbar;
c.Label.String = 'Delay (ms)';
xlabel('Electrodes');
ylabel('Electrodes');

% load electrode pattern
load('./computed_sphere_grid/node_pos_19.mat');

% find center electrode
[~, center] = min((position(:,1) - mean(position(:,1))).^2 + (position(:,2) - mean(position(:,2))).^2);

figure;
scatter(position(:,1), position(:,2), 400, 1000 * delay(center,:), 'filled');
hold on
plot(position(center,1), position(center,2), 'X');
c = colorbar;
c.Label.String = 'Delay to center electrode X (ms)';
xlabel('Electrode position (mm)');
ylabel('Electrodes position (mm)');


%% Estimate parameters of the waves
[src_dir, speed, ci_dir, ci_sp] = estimate_wave(delay, position, 'plot');


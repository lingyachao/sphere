%% if drawing one period, draw this
central_t = 5:10:195;
period_idx = find(central_t == 115);

%% coherence between two farthest electrodes from 110-120s
figure;
subplot(2, 2, 1);
plot(freq, squeeze(t_coh(10,13,:,period_idx)));
xlabel('frequency (Hz)');
ylabel('coherence');
title('coherence between two farthest electrodes during 110-120s');

%% average coherence through time

% load electrode positions
load('./computed_sphere_grid/node_pos_13.mat');
% position = locs(macro_idx(locs(macro_idx,1) > 45), [3,2]);

ic = NaN(P,1);
for i = 1:P
    coh_avg = mean(squeeze(t_coh(:,:,:,i)), 3);
    dist = pdist(position)';
    coh_pair = tril(coh_avg, -1);
    coh_pair = coh_pair(:);
    coh_pair = coh_pair(coh_pair > 0);
    f = polyfit(dist, coh_pair, 1);
    ic(i) = f(2);
end

subplot(2, 2, 2);
plot(central_t, ic);
line([30 30], [-100 100], 'Color', 'r', 'LineWidth', 2);
line([150 150], [-100 100], 'Color', 'r', 'LineWidth', 2);
xlabel('time (s)');
xlim([0 200]);
ylabel('average coherence (intercept at 0mm)');
ylim([0 1.2]);
title('average coherence through time');

%% estimate delays between electrodes during 110-120s

[delay, delay_ci_lo, delay_ci_up] = compute_delay(t_coh(:,:,:,period_idx), ...
    t_coh_conf(:,:,:,period_idx), t_phi(:,:,:,period_idx), freq);

% find center electrode
[~, center] = min((position(:,1) - mean(position(:,1))).^2 + (position(:,2) - mean(position(:,2))).^2);

% subplot(2, 2, 3);
% scatter(position(:,1), position(:,2), 400, 1000 * delay(center,:), 'filled');
% c = colorbar;
% c.Label.String = 'delay to center electrode X (ms)';
% xlabel('electrode position (mm)');
% ylabel('electrodes position (mm)');
% title('delays to central electrode during 110-120s');

%% estimate parameters of the wave during 110-120s

subplot(2, 2, 3);
[src_dir, speed, ci_dir, ci_sp] = estimate_wave(delay, position, 'plot');
c = colorbar;
c.Label.String = 'delay to center electrode X (ms)';
title('fit a plane wave (during 110-120s)');

%% estimate directions of all periods

dirs = NaN(P,1);
for i = 1:P
    [delay, delay_ci_lo, delay_ci_up] = compute_delay(t_coh(:,:,:,i), ...
        t_coh_conf(:,:,:,i), t_phi(:,:,:,i), freq);
    [src_dir, speed, ci_dir, ci_sp] = estimate_wave(delay, position, 'no');
    dirs(i) = src_dir;
end

subplot(2, 2, 4);
scatter(central_t, dirs, 10, 'filled');
line([30 30], [-100 100], 'Color', 'r', 'LineWidth', 2);
line([150 150], [-100 100], 'Color', 'r', 'LineWidth', 2);
xlabel('time (s)');
xlim([0 200]);
ylabel('source direction (rad)');
ylim([-pi pi]);
title('source direction through time');
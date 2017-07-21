%% coherence between two farthest electrodes from 110-120s

fg = figure;
subplot(2, 2, 1);
plot(freq, squeeze(t_coh(1,8,:,period_idx)));
xlabel('frequency (Hz)');
ylabel('coherence');
title('coherence between two farthest electrodes during 110-120s');

%% average coherence through time

ic = NaN(P,1);
for i = 1:P
    coh_avg = mean(squeeze(t_coh(:,:,:,i)), 3);
    dist_pairs = pdist(electrode_2d)';
    coh_pair = tril(coh_avg, -1);
    coh_pair = coh_pair(:);
    coh_pair = coh_pair(coh_pair > 0);
    f = polyfit(dist_pairs, coh_pair, 1);
    ic(i) = f(2);
end

subplot(2, 2, 2);
plot(central_t, ic);
xlabel('time (s)');
xlim([0 200]);
ylabel('average coherence (intercept at 0mm)');
ylim([0 1.2]);
title('average coherence through time');

%% estimate delays between electrodes during 110-120s

[delay, delay_ci_lo, delay_ci_up] = compute_delay(t_coh(:,:,:,period_idx), ...
    t_coh_conf(:,:,:,period_idx), t_phi(:,:,:,period_idx), freq);

% find center electrode
[~, center] = min((electrode_2d(:,1) - mean(electrode_2d(:,1))).^2 + (electrode_2d(:,2) - mean(macro_2d(:,2))).^2);

%% estimate parameters of the wave during 110-120s

subplot(2, 2, 3);
[src_dir, speed, ci_dir, ci_sp] = estimate_wave(delay, electrode_2d, 'plot');
c = colorbar;
c.Label.String = 'delay to center electrode X (ms)';
title('fit a plane wave (during 110-120s)');

%% estimate directions of all periods

dirs = NaN(P,1);
speeds = NaN(P,1);
for i = 1:P
    [delay, ~, ~] = compute_delay(t_coh(:,:,:,i), ...
        t_coh_conf(:,:,:,i), t_phi(:,:,:,i), freq);
    [dirs(i), speeds(i), ~, ~] = estimate_wave(delay, electrode_2d, 'no');
end

subplot(2, 2, 4);
scatter(central_t, dirs, 10, 'filled');
ylabel('source direction (rad)');
ylim([-pi pi]);

yyaxis right;
scatter(central_t, speeds, 10, 'd');
ylabel('velocity (cm/s)');

xlabel('time (s)');
xlim([0 200]);
title('source direction through time');

%% save figure
saveas(fg, COHERENCE_FIG);

%% display average speed in the last 5 periods
fprintf(['SWD speed is ' num2str(mean(speeds(end-4:end))) ' cm/s \n']);

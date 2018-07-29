%% coherence between two farthest electrodes chosen period
subplot(2, 2, 1);
plot(freq, squeeze(t_coh(1,8,:,period_idx)));
xlabel('frequency (Hz)');
ylabel('coherence');

period_start = central_t(period_idx) - wind/2;
period_str = [num2str(period_start) '-' num2str(period_start + wind) 's'];
% title(['coherence between two farthest electrodes during ' period_str]);

%% average coherence through time
ic = NaN(P,1);
for i = 1:P
    coh_avg = mean(squeeze(t_coh(:,:,:,i)), 3);
    dist_pairs = pdist(electrode_2d)';
    coh_pair = tril(coh_avg, -1);
    coh_pair = coh_pair(:);
    coh_pair = coh_pair(coh_pair > 0 | isnan(coh_pair));
    f = polyfit(dist_pairs, coh_pair, 1);
    ic(i) = f(2);
end

subplot(2, 2, 4);
plot(central_t, ic);
xlabel('time (s)');
xlim([25 total_time]);
ylabel('global coherence');
ylim([0 1.2]);
% title('average coherence through time');

%% estimate delays between electrodes during chosen period
[delay, delay_ci_lo, delay_ci_up] = compute_delay(t_coh(:,:,:,period_idx), ...
    t_coh_conf(:,:,:,period_idx), t_phi(:,:,:,period_idx), freq);

% find center electrode
[~, center] = min((electrode_2d(:,1) - mean(electrode_2d(:,1))).^2 + (electrode_2d(:,2) - mean(macro_2d(:,2))).^2);

%% estimate parameters of the wave during chosen period
sp = subplot(2, 2, 2);
[src_dir, speed, ci_dir, ci_sp] = estimate_wave(delay, electrode_2d, 'plot');
c = colorbar(sp);
c.Label.String = 'delay to center electrode X (ms)';
% title(['fit a plane wave (during ' period_str ')']);

%% estimate directions of all periods
dirs = NaN(P,1);
speeds = NaN(P,1);
for i = 1:P
    [delay, ~, ~] = compute_delay(t_coh(:,:,:,i), ...
        t_coh_conf(:,:,:,i), t_phi(:,:,:,i), freq);
    [dirs(i), speeds(i), ~, ~] = estimate_wave(delay, electrode_2d, 'no');
end

subplot(2, 2, 3);
scatter(central_t, dirs, 10, 'filled');
ylabel('source direction (rad)');
ylim([-pi pi]);
xlim([25 total_time]);

yyaxis right;
scatter(central_t, speeds, 10, 'd');
ylabel('velocity (cm/s)');

xlabel('time (s)');
xlim([25 total_time]);
% title('source direction through time');

%% display average speed in the last 3 periods
swd_speed = mean(speeds(end-7:end-5));
fprintf(['SWD speed is ' num2str(swd_speed) ' cm/s \n']);

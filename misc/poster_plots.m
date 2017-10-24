figure;
set(0, 'DefaultAxesFontSize', 8);

%% single node
subplot(3, 1, 1);

plot(sparse_time, single_node(:,[1,2,9]));
xlabel('time (s)');
xlim([30 total_time]);
ylabel('firing rate (Hz)');
legend('E', 'I', 'I\_{fs}');
title('firing rate of each population vs. time');


%% coherence
subplot(3, 1, 2);

central_t = int32(total_time * (1/P/2 : 1/P : 1-1/P/2));
    period_idx = length(central_t) - 2;
[t_coh,t_coh_conf,t_phi,electrode_2d] = deal( ...
    macro_t_coh, macro_t_coh_conf, macro_t_phi, macro_2d);

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

plot(central_t, ic);
xlabel('time (s)');
xlim([30 total_time]);
ylabel('average coherence');
ylim([0 1.2]);
title('average coherence vs. time');

%% direction

[t_coh,t_coh_conf,t_phi,electrode_2d] = deal( ...
    micro_t_coh, micro_t_coh_conf, micro_t_phi, micro_2d);

[delay, delay_ci_lo, delay_ci_up] = compute_delay(t_coh(:,:,:,period_idx), ...
    t_coh_conf(:,:,:,period_idx), t_phi(:,:,:,period_idx), freq);

% find center electrode
[~, center] = min((electrode_2d(:,1) - mean(electrode_2d(:,1))).^2 + (electrode_2d(:,2) - mean(macro_2d(:,2))).^2);

% estimate directions of all periods
dirs = NaN(P,1);
speeds = NaN(P,1);
for i = 1:P
    [delay, ~, ~] = compute_delay(t_coh(:,:,:,i), ...
        t_coh_conf(:,:,:,i), t_phi(:,:,:,i), freq);
    [dirs(i), speeds(i), ~, ~] = estimate_wave(delay, electrode_2d, 'no');
end

subplot(3, 1, 3);
scatter(central_t, dirs, 10, 'filled');
ylabel('source direction (rad)');
ylim([-pi pi]);

yyaxis right;
scatter(central_t, speeds, 10, 'd');
ylabel('velocity (cm/s)');

xlabel('time (s)');
xlim([30 total_time]);
title('source direction / traveling wave velocity vs. time');
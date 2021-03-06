clear; close all;
clearvars -global;

%% specify run type
DATA_STORAGE = 'C:/Users/monica/simulation_data/';
type = 'sphere';
note = 'Nie670';

save_output = false;
visualize = true;
print_count = true;
use_fluc = false;
flag_dense_coh = false;

K = ~use_fluc * 4000 + use_fluc * 10000;

%% run setup script
setup;

%% make modifications to constants

% HL.Nie_fs = 10;
% HL.theta_e = -58.8;
% HL.sigma_e = 3.0;
% HL.Qe_max = 35;

if strcmp(type, 'sphere')
    HL.kR = 5;
    HL.k_decay = 0.01;
    HL.KtoD = -1.5;
else
    HL.kR = 10;
    HL.k_decay = 1;
    HL.KtoD  = -2.5;
end

HL.k_decay = HL.k_decay * ones(N,1);
% HL.k_decay = HL.k_decay * max((1 + 0.5*randn(N,1)), 0);
% HL.k_decay = HL.k_decay / 2 + (laplacian .* (laplacian > 0)) * HL.k_decay / 12;
% HL.k_decay(zones.focus_zone) = 0.001;

% low_kdecay_idx = randsample(find(zones.lessihb_zone), 20);
% HL.k_decay(sum(abs(laplacian(:,low_kdecay_idx)), 2) > 0) = 0.01;
HL.k_decay(zones.normal_zone) = 100;

HL.prodRatio = 0.1;
HL.KtoVe = 0; HL.KtoVi = 0; HL.KtoVi_fs = 0;
HL.D22min = 0.1;

last.D22(:) = 5; last.D11 = last.D22/100;
last.K(:) = 5;
if use_fluc
    last.K(zones.focus_zone | zones.lessihb_zone) = 5.5;
end

% seizure initiation source
if ~use_fluc
    source_drive_amp = 8;
    source_drive_duration = 10; % seconds
else
    % distributions for fluctuations
    sc_basal = 1;
    sc_rand_width = 1;
    sc_high_diff_max = 28;
    sc_high_prob = 0.0005;
    pd = makedist('Binomial', 'N', 1, 'p', sc_high_prob);
    phi_ee_sc_base = HL.phi_ee_sc(1);
end

if save_output
    save(META_FILE, 'HL', 'map', 'fine_idx', 'single_node_idx', 'last');
end

[trace.t, trace.Qe, trace.Qi, trace.Qi_fs, trace.K, trace.Ve] = deal(NaN);

%% run simulation
for k = 1:K
     
    if ~use_fluc && k < source_drive_duration/T0
        source_drive = source_drive_amp;
    else
        source_drive = NaN;
    end

    if print_count
        fprintf(['Running simulation , ' num2str(k) ' ... ']);
    end
    tic;
  
    if use_fluc && rem(k, 2) == 1
        
        % sc_high_diff = sc_high_diff_max;
        sc_high_diff = sc_high_diff_max * exp(-4 * mean(last.Qe));
        
        HL.phi_ee_sc = sc_basal + sc_rand_width * randn(N, 1) + ...
            sc_high_diff * random(pd, N, 1);
        HL.phi_ee_sc = max(HL.phi_ee_sc, 0);

        HL.phi_ee_sc = HL.phi_ee_sc * phi_ee_sc_base;
    end
    
    [samp_time,last,fine] = seizing_cortical_field(...
        source_drive, map, T0, last, ...
        locs, laplacian, laplacian, avg_D, ...
        zones, fine_idx, single_node_idx, ...
        save_output);
    
    if visualize
        clf(f);
        if strcmp(type, 'sphere')
            trace.t = [trace.t k*T0];
            trace.Qe = [trace.Qe last.Qe(100)];
            trace.Qi = [trace.Qi last.Qi(100)];
            trace.Qi_fs = [trace.Qi_fs last.Qi_fs(100)];
            trace.K = [trace.K last.K(100)];
            trace.Ve = [trace.Ve last.Ve(100)];
            plot_sphere_instance(locs, last, NaN, NaN, trace);
        elseif strcmp(type, 'brain')
            plot_brain_instance(surf, surf_sphere, last, NaN, NaN);
        end
        drawnow;
    end
    
    if save_output
        save([OUTPUT_DIR 'seizing_cortical_field_k_' num2str(k) '.mat'], ...
            'samp_time', 'last', 'fine');
    end

    if print_count
        fprintf(['RT ' num2str(toc) '\n']);
    end
end

%% run analysis
if save_output
    main_plot_graphs(id, DATA_STORAGE, true, false, true, flag_dense_coh);
end
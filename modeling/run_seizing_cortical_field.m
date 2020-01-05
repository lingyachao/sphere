clear; close all;
clearvars -global;

%% specify run type
DATA_STORAGE = 'C:/Users/monica/simulation_data/';
type = 'sphere';
note = 'Nie670_kdecay0.1';

save_output = true;
visualize = true;
print_count = true;
use_fluc = false;
flag_dense_coh = true;

K = ~use_fluc * 4000 + use_fluc * 10000;

%% run setup script
setup;

%% make modifications to constants
if strcmp(type, 'sphere')
    HL.kR = 5;
    HL.k_decay = 0.1;
    HL.KtoD = -1.5;
else
    HL.kR = 10;
    HL.k_decay = 1;
    HL.KtoD  = -2.5;
end

HL.k_decay = HL.k_decay * ones(N,1);
HL.k_decay(zones.normal_zone) = 100;

HL.prodRatio = 0.1;
HL.KtoVe = 0; HL.KtoVi = 0; HL.KtoVi_fs = 0;
HL.D22min = 0.1;

last.D22(:) = 5; last.D11 = last.D22/100;
last.K(:) = 5;
if use_fluc
    last.K(zones.focus_zone | zones.lessihb_zone) = 5.4;
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
    sc_high_prob = 0.001;
    pd = makedist('Binomial', 'N', 1, 'p', sc_high_prob);
    phi_ee_sc_base = HL.phi_ee_sc(1);
end

if save_output
    save(META_FILE, 'HL', 'map', 'fine_idx', 'single_node_idx', 'last');
end

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
        
        sc_high_diff = sc_high_diff_max;
        % sc_high_diff = sc_high_diff_max * exp(-3 * mean(last.Qe));
        
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
            plot_sphere_instance(locs, last, NaN, NaN);
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
        % fprintf(['RT ' num2str(toc) ' sc_high_diff ' num2str(sc_high_diff) '\n']);
    end
end

%% run analysis
if save_output
    main_plot_graphs(id, DATA_STORAGE, true, false, true, flag_dense_coh);
end
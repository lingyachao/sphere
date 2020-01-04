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

K = ~use_fluc * 3000 + use_fluc * 20000;

%% load grid
if strcmp(type, 'sphere')
    load('N10242_R10.mat'); 
elseif strcmp(type, 'brain')
    load('N40962.mat'); avg_D = 0.3777;
    load('unitsphere.mat');
    
    surf.vertices = locs;
    surf.faces = tri;
    surf_sphere.vertices = 10 * coord';
    surf_sphere.faces = tri;
else
    error('not recognized type');
end

%% set plotting window
if visualize
    f = figure;
    set(f, 'Position', [200 300 900 400]);
end

%% initialize initial state
last = make_IC(N);
last.Qi_fs = last.Qi;
last.Vi_fs = last.Vi;
last.F_ii_fs = last.F_ii;
last.Phi_ii_fs = last.Phi_ii;
last.dVi_fs = last.dVi;

% [DATA_DIR, type] = find_full_id('./data/', '09211929');
% load([DATA_DIR 'raw/seizing_cortical_field_k_1500.mat'], 'last');

%% define zones
map = make_map(laplacian);

if strcmp(type, 'sphere')
    lessihb_filter = lessihb_area;
else
    lessihb_filter = lessihb_area_brain;
end

zones.focus_zone = map == 1;
zones.lessihb_zone = lessihb_filter & map ~= 1;
zones.normal_zone = ~lessihb_filter;

if strcmp(type, 'sphere')
    macro_idx = [744 437 821 1141 1140 820 436 251 555 981 1253 1585 1537 1584 1252 980 554 250 187];
    micro_idx = [744 659 753 837 836 752 658 579 669 777 845 933 929 932 844 776 668 578 573];
    fine_idx = union(find(map), [macro_idx, micro_idx]);
else
    fine_idx = union(find(map), find(coord(1,:)' > 0.7));
end

% single node index
if strcmp(type, 'sphere')
    single_node_idx = 744;
else
    single_node_idx = 29584;
end
 
%% initialize constants and make modifications
T0 = 0.1;

global HL
HL = SCM_init_globs(N);

if strcmp(type, 'sphere')
    HL.kR = 5;
    HL.k_decay = 0.1;
    HL.KtoD = -1.5;    
else
    HL.kR = 10;
    HL.k_decay = 1;
    HL.KtoD  = -2.5;
end

% HL.kR = HL.kR * ones(N, 1);
% HL.kR(zones.normal_zone) = 0;

HL.k_decay = HL.k_decay * ones(N,1);
HL.k_decay(zones.normal_zone) = 100;

HL.prodRatio = 0.1;
HL.KtoVe = 0;
HL.KtoVi = 0;
HL.KtoVi_fs = 0;
HL.D22min = 0.1;
HL.FS_ratio = 0;

last.D22(:) = 5; last.D11 = last.D22/100;
last.K(:) = 5;
if use_fluc
    last.K(zones.focus_zone | zones.lessihb_zone) = 5.4;
end

if ~use_fluc
    % last.K(map == 1) = 12;
    source_drive_amp = 8;
    source_drive_duration = 10; % seconds
else
    % HL.Nie_fs = HL.Nie_fs * ones(N, 1);
    % HL.Nie_fs(map == 1) = 0;
    phi_ee_sc_base = HL.phi_ee_sc(1);
end

% distributions for fluctuations
sc_basal = 1;
sc_rand_width = 1;
sc_high_diff = 28;
sc_high_prob = 0.001;
pd = makedist('Binomial', 'N', 1, 'p', sc_high_prob);

%% set the output directory and save meta file
if save_output
    id = datestr(now, 'yy_mmdd_HHMM');
    if strcmp(type, 'sphere')
        folder_name = ['sphere_N' num2str(N) '_R' num2str(R) '_' id '_' note];
    elseif strcmp(type, 'brain')
        folder_name = ['brain_N' num2str(N) '_' id '_' note];
    end

    OUTPUT_DIR = [DATA_STORAGE folder_name '/raw/'];
    mkdir(OUTPUT_DIR);
    
    META_FILE = [DATA_STORAGE folder_name '/vars.mat'];
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
        
        HL.phi_ee_sc = sc_basal + sc_rand_width * randn(N, 1) + ...
            sc_high_diff * random(pd, N, 1);
        HL.phi_ee_sc = max(HL.phi_ee_sc, 0);

        % if all(last.Qe(map == 1) < 3 & last.K(map == 1) < 7)
        %    HL.phi_ee_sc(map == 1) = datasample(tail_discrete, 7, 'Weight', pdfs);
        %    HL.phi_ee_sc(map == 1)'
        % end
        
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
    end
end

%% run analysis
if save_output
    main_plot_graphs(id, DATA_STORAGE, true, false, true, flag_dense_coh);
end
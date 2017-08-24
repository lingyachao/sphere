clear; close all;

%% specify run type
type = 'sphere';
note = 'depolarization_realrest_maxK12_onlyFSproduceK_nodecay';
save_output = true;
visualize = true;
print_count = true;

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

%% initialize parameters and map
K = 3000;
T0 = 0.1;
map = make_map(laplacian);

%% initialize initial state
last = make_IC(N);
last.Qi_fs = last.Qi;
last.Vi_fs = last.Vi;
last.F_ii_fs = last.F_ii;
last.Phi_ii_fs = last.Phi_ii;
last.dVi_fs = last.dVi;
% load('./data/sphere_N10242_R10_08012157_depolarization_2pops_activation20/raw/seizing_cortical_field_k_900.mat', 'last');

%% define zones
lessihb_filter = true(N, 1);

zones.focus_zone = map == 1;
zones.lessihb_zone = lessihb_filter & map ~= 1;
zones.normal_zone = ~lessihb_filter;

% indices of normal-zone nodes to be recorded at fine time-scale
normal_sample_idx = []; % randsample(find(zones.normal_zone),3);

if strcmp(type, 'sphere')
    macro_idx = [744 437 821 1141 1140 820 436 251 555 981 1253 1585 1537 1584 1252 980 554 250 187];
    micro_idx = [744 659 753 837 836 752 658 579 669 777 845 933 929 932 844 776 668 578 573];
    fine_idx = union(find(map), [macro_idx, micro_idx]);
else
    fine_idx = union(find(map), find(coord(1,:)' > 0.95));
end
 
%% initialize constants and make modifications
global HL
HL = SCM_init_globs(N);

HL.kR = 15 * ones(N,1);
% HL.kR(zones.normal_zone) = 0;

HL.k_decay = 0;

HL.KtoVe = 0;
HL.KtoVi = 0;
HL.KtoVi_fs = 0;
HL.KtoD  = -5;
HL.D22min = 0.1;
HL.FS_ratio = 0;

% [HL.Nee_a, HL.Nei_a] = deal(1000, 1000);
last.D22(:) = 7; last.D11 = last.D22/100;
last.K(:) = 5;
% last.dVe(:) = -3;
% last.dVi(:) = 0;

%% set the output directory and save meta file
if save_output
    id = datestr(now, 'mmddHHMM');
    if strcmp(type, 'sphere')
        folder_name = ['sphere_N' num2str(N) '_R' num2str(R) '_' id '_' note];
    elseif strcmp(type, 'brain')
        folder_name = ['brain_N' num2str(N) '_' id '_' note];
    end

    OUTPUT_DIR = ['./data/' folder_name '/raw/'];
    mkdir(OUTPUT_DIR);
    
    META_FILE = ['./data/' folder_name '/vars.mat'];
    save(META_FILE, 'HL', 'map', 'fine_idx', 'normal_sample_idx', 'last');
end

%% run simulation
for k = 1:K
     
    if true % k < 2000
        source_drive = 5;
    else
        source_drive = NaN;
    end

    if print_count
        fprintf(['Running simulation , ' num2str(k) ' ... ']);
    end
    tic;

    [samp_time,last,fine] = seizing_cortical_field(...
        source_drive, map, T0, last, ...
        locs, laplacian, avg_D, ...
        zones, fine_idx, normal_sample_idx, ...
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
        % fprintf(['K at node8 ' num2str(last.K(8)) '\n']);
        % fprintf(['dVi at node8 ' num2str(last.dVi(8)) '\n']);
        % fprintf(['Vi at node8 ' num2str(last.Vi(8)) '\n']);
    end
end

%% run analysis
if save_output
    main_plot_graphs(id, './data/', true, false, true);
end
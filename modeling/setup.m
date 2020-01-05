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

%% initialize constants
global HL
HL = SCM_init_globs(N);
T0 = 0.1;

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

%% set the output directory
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
end
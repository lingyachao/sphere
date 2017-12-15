clear; close all;

%% specify run type
type = 'sphere';
note = 'depolarization_realrest_maxK12_withnormal_centerK12_randomarea_changerev_changeNa_Nie665';
save_output = true;
visualize = true;
print_count = true;

%% load grid
if strcmp(type, 'sphere')
    load('N10242_R10.mat'); 
    % avg_D = avg_D + 0.5 * (rand(N, 1));
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

% [DATA_DIR, type] = find_full_id('./data/', '09211929');
% load([DATA_DIR 'raw/seizing_cortical_field_k_1500.mat'], 'last');

%% define zones
if strcmp(type, 'sphere')
    lessihb_filter = lessihb_area;
else
    lessihb_filter = lessihb_area_brain;
end

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
    fine_idx = union(find(map), find(coord(1,:)' > 0.7));
end
 
%% initialize constants and make modifications
global HL
HL = SCM_init_globs(N);

if strcmp(type, 'sphere')
    HL.kR = 5;
    HL.k_decay = 0;
    HL.KtoD = -1.5;    
else
    HL.kR = 10;
    HL.k_decay = 1;
    HL.KtoD  = -2.5;
end

HL.kR = HL.kR * ones(N, 1);
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
last.K(map == 1) = 12;

% HL.Nie_fs = HL.Nie_fs * ones(N, 1);
% HL.Nie_fs(map == 1) = 0;
phi_ee_sc_base = HL.phi_ee_sc(1);

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

%% tail of normal distribution
tail_discrete = (18:0.1:50)';
pdfs = normpdf(tail_discrete, 1, 1);

%% run simulation
for k = 1:K
     
    if true % k < 2000
        source_drive = NaN;
    else
        source_drive = NaN;
    end

    if print_count
        fprintf(['Running simulation , ' num2str(k) ' ... ']);
    end
    tic;
     
%     if k < 5
%         HL.phi_ee_sc(laplacian(:,200) ~= 0) = 30 * phi_ee_sc_base;
%     else
%         HL.phi_ee_sc(laplacian(:,200) ~= 0) = phi_ee_sc_base;
%     end
    
%     HL.phi_ee_sc = randn(N, 1)*1 + 1;
%     HL.phi_ee_sc = max(HL.phi_ee_sc, 0);
%     
%     if all(last.Qe(map == 1) < 3 & last.K(map == 1) < 7)
%         HL.phi_ee_sc(map == 1) = datasample(tail_discrete, 7, 'Weight', pdfs);
%         HL.phi_ee_sc(map == 1)'
%     end
    
%     HL.phi_ee_sc = HL.phi_ee_sc * phi_ee_sc_base;
    
    [samp_time,last,fine] = seizing_cortical_field(...
        source_drive, map, T0, last, ...
        locs, laplacian, laplacian, avg_D, ...
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
        fprintf(['K at node2 ' num2str(last.K(2)) '\n']);
        fprintf(['K at node500 ' num2str(last.K(500)) '\n']);
        % fprintf(['dVi at node8 ' num2str(last.dVi(8)) '\n']);
        % fprintf(['Vi at node8 ' num2str(last.Vi(8)) '\n']);
    end
end

%% run analysis
if save_output
    main_plot_graphs(id, './data/', true, false, true);
end
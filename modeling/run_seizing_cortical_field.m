clear; close all;

%% specify run type
type = 'sphere';
note = 'full_FSinhib';
save_output = true;
visualize = true;
print_count = true;

%% load grid
if strcmp(type, 'sphere')
    load('N10242_R10.mat'); 
elseif strcmp(type, 'brain')
    load('N40962.mat');
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
k = 0;
K = 2000;
T0 = 0.1;
map = make_map(laplacian);

%% initialize initial state
last = make_IC(N);

%% define zones
lessihb_filter = true(N, 1);

zones.focus_zone = map == 1;
zones.lessihb_zone = lessihb_filter & map ~= 1;
zones.normal_zone = ~lessihb_filter;

lessihb_idx = find(lessihb_filter);
normal_sample_idx = []; % randsample(find(zones.normal_zone),3);
    
%% initialize constants and make modifications
global HL
HL = SCM_init_globs(N);

HL.kR = 2.5;
HL.KtoVe = 0;
HL.KtoVi = 0;
HL.KtoD  = -20;
HL.D22min = 0.1;

% last.D22(:) = 0.5;
% last.D11 = last.D22/100;
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
    save(META_FILE, 'HL', 'map', 'lessihb_idx', 'normal_sample_idx', 'last');
end

%% run simulation
for k = 1:K
     
    if true
        source_drive = 3;
    elseif k > 150 / T0
        source_drive = NaN;
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
        zones, lessihb_idx, normal_sample_idx, ...
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
        % fprintf(['mean ' num2str(mean(last.Ve)) ' sd ' num2str(std(last.Ve)) '\n']);
        % fprintf(['K normal ' num2str(mean(last.K(zones.normal_zone))) ' K abnormal ' num2str(mean(last.K(lessihb_idx))) '\n']);
        % fprintf(['D2 ' num2str(mean(last.D22(lessihb_idx))) ' dVe ' num2str(mean(last.dVe(lessihb_idx))) '\n']);
        % fprintf(['Ve focus ' num2str(last.Ve(1)) '\n']);
    end
end

%% run analysis
if save_output
    main_plot_graphs;
end
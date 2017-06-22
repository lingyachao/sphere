clear; close all;

%% specify run type
type = 'sphere';
note = 'full_data';
save_output = false;
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
% lessihb_filter = lessihb_area;
lessihb_filter = locs(:,3) < -6;
% lessihb_filter = coord(1,:)' > 0.5;
% lessihb_filter = true(N, 1);

zones.focus_zone = map == 1;
zones.lessihb_zone = lessihb_filter & map ~= 1;
zones.normal_zone = ~lessihb_filter;

lessihb_idx = find(lessihb_filter);
normal_sample_idx = randsample(find(zones.normal_zone),3);
    
%% initialize constants and make modifications
global HL
HL = SCM_init_globs(N);

% HL.kR = 0;
HL.KtoD  = 0;

last.D22(:) = 2;
last.D11 = last.D22/100;
% last.dVe(:) = -8;
% last.dVi(:) = 0;

% HL.FSi = [];
HL.FSi = false(N, 1);
HL.FSi(randsample(1:N, floor(N/5))) = true;

last.Qi_fs = last.Qi;
last.Qi_fs(~HL.FSi) = 0;

% increase inhibitory strength in all locations other than a patch
HL.Vi_rest(zones.normal_zone) = HL.Vi_rest(zones.normal_zone) + 3;

%% 
% Nie_b(normal_zone) = 1.2 * HL.Nie_b;
% Nii_b(normal_zone) = 1.2 * HL.Nii_b;
% Nie_b(focus_zone) = 0.95 * HL.Nie_b;
% Nii_b(focus_zone) = 0.95 * HL.Nii_b;
% Nie_b(lessihb_zone) = 0.95 * HL.Nie_b;
% Nii_b(lessihb_zone) = 0.95 * HL.Nii_b;

% [HL.Nee_a, HL.Nei_a] = deal(5000, 5000);
% lam = 1;
% HL.gamma_i = HL.gamma_i / lam;
% HL.gi = HL.gi * lam;

%% set the output directory and save meta file
if save_output
    if strcmp(type, 'sphere')
        folder_name = ['sphere_N' num2str(N) '_R' num2str(R) '_' datestr(now, 'mmddHHMM') '_' note];
    elseif strcmp(type, 'brain')
        folder_name = ['brain_N' num2str(N) '_' datestr(now, 'mmddHHMM') '_' note];
    end

    OUTPUT_DIR = ['./data/' folder_name '/raw/'];
    mkdir(OUTPUT_DIR);
    
    META_FILE = ['./data/' folder_name '/meta.mat'];
    save(META_FILE, 'HL', 'map', 'lessihb_idx', 'normal_sample_idx', 'last');
end

%% run simulation
for k = 1:K
     
    if true
        source_drive = 2;
    elseif k > 150 / T0
        source_drive = NaN;
    else
        source_drive = NaN;
    end
    
    % HL.Qi_max = HL.Qi_max - 0.2;

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

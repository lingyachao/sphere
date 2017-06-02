clear; close all;

%% specify run type
type = 'sphere';
save_output = false;

%% load grid
if strcmp(type, 'sphere')
    % load('./computed_sphere_grid/N10242_R10_wideNodes.mat');
    load('./computed_sphere_grid/N42_R10.mat');
    pos_hemi = locs(:,3) >= 0;
    neg_hemi = locs(:,3) < 0;
elseif strcmp(type, 'brain')
    load('./computed_brain_grid/N40962.mat');
    load('unitsphere.mat', 'coord');
else
    error('not recognized type');
end

%% set the output directory
note = 'no_stimuli_noK';

if save_output
    if strcmp(type, 'sphere')
        folder_name = ['sphere_N' num2str(N) '_R' num2str(R) '_' datestr(now, 'mmddHHMM') '_' note];
    elseif strcmp(type, 'brain')
        folder_name = ['brain_N' num2str(N) '_' datestr(now, 'mmddHHMM') '_' note];
    end

    OUTPUT_DIR = ['./data/' folder_name '/raw/'];
    mkdir(OUTPUT_DIR);
end

%% set plotting window
f = figure;
set(f, 'Position', [200 300 900 400]);
if strcmp(type, 'brain')
    load('autism.surface.mat', 'tri');
    surf.vertices = locs;
    surf.faces = tri;
end

%% initialize parameters and map
k = 0;
K = 10;
T0 = 0.1;
visualize = true;
map = make_map(laplacian);

%% initialize initial state
last = make_IC(N);
last.dVe = zeros(N, 1);
last.dVi = zeros(N, 1);
% last.Ve = -60 * ones(N, 1);

%% define zones
% lessihb_idx = lessihb_area;
% lessihb_idx = locs(:,3) < -6;
% lessihb_idx = coord(1,:)' > 0.5;
lessihb_idx = false(N, 1);
zones.focus_zone = map == 1;
zones.lessihb_zone = lessihb_idx & map ~= 1;
zones.normal_zone = ~lessihb_idx;
    
%% initialize constants and make modifications
global HL
HL = SCM_init_globs;

% increase inhibitory strength in all locations other than a patch
HL.Nie_b = HL.Nie_b * ones(N, 1);
HL.Nii_b = HL.Nii_b * ones(N, 1);
HL.Vi_rest = HL.Vi_rest * ones(N, 1);
       
% Nie_b(normal_zone) = 1.2 * HL.Nie_b;
% Nii_b(normal_zone) = 1.2 * HL.Nii_b;
% Nie_b(focus_zone) = 0.95 * HL.Nie_b;
% Nii_b(focus_zone) = 0.95 * HL.Nii_b;
% Nie_b(lessihb_zone) = 0.95 * HL.Nie_b;
% Nii_b(lessihb_zone) = 0.95 * HL.Nii_b;

% Vi_rest(normal_zone) = Vi_rest(normal_zone) + 3;

% [HL.Nee_a, HL.Nei_a] = deal(5000, 5000);
% lam = 1;
% HL.gamma_i = HL.gamma_i / lam;
% HL.gi = HL.gi * lam;

%% run simulation
for k = 1:K
    
    if k <= 30 / T0 % k == 20 / T0
        source_drive = NaN;
    elseif k > 150 / T0
        source_drive = NaN;
    else
        source_drive = NaN; % 3;
    end

    fprintf(['Running simulation , ' num2str(k) ' ... ']);
    tic; % start timer

    [samp_time,last,fine] = seizing_cortical_field(...
        source_drive, map, T0, last, ...
        locs, laplacian, avg_D, ...
        zones, micro_idx, macro_idx, focus_idx, normal_idx);
    
    fprintf(['RT ' num2str(toc) '\n']);

    if visualize
        
        if strcmp(type, 'sphere')
            subplot(1, 2, 1);
            scatter(locs(pos_hemi,1), locs(pos_hemi,2), 15, last.Qe(pos_hemi), 'filled');
            caxis([0,30]);

            subplot(1, 2, 2);    
            scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
            caxis([0,30]);
            
        elseif strcmp(type, 'brain')
            clf(f);
            figure_wire(surf, last.Qe, false);
            
        end
        
        drawnow;
    end
    
    if save_output
        save([OUTPUT_DIR 'seizing_cortical_field_k_' num2str(k) '.mat'], ...
            'samp_time', 'last', 'fine');
    end
    
    fprintf(['mean ' num2str(mean(last.Ve)) ' sd ' num2str(std(last.Ve)) '\n']);
end 

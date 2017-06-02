clear; close all;

type = 'sphere';
save_output = false;

if strcmp(type, 'sphere') % load sphere grid
    % load('./computed_sphere_grid/N10242_R10_wideNodes.mat');
    load('./computed_sphere_grid/N42_R10.mat');
    pos_hemi = locs(:,3) >= 0;
    neg_hemi = locs(:,3) < 0;
elseif strcmp(type, 'brain') % load brain grid
    load('./computed_brain_grid/N40962.mat');
    load('unitsphere.mat', 'coord');
else
    error('not recognized type');
end

% set the output directory
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

% initialize parameters, initial state and map
k = 0;
K = 200;
T0 = 0.1;
visualize = true;
map = make_map(laplacian);
last = make_IC(N);

last.dVe = zeros(N, 1); %ones(N, 1) * 2;
last.dVi = zeros(N, 1);
last.Ve = -60 * ones(N, 1);

% define nodes that have less inhibition
% lessihb_idx = lessihb_area;
% lessihb_idx = locs(:,3) < -6;
% lessihb_idx = coord(1,:)' > 0.5;
lessihb_idx = false(N, 1);

% set plotting window
f = figure;
set(f, 'Position', [200 300 900 400]);
if strcmp(type, 'brain')
    load('autism.surface.mat', 'tri');
    surf.vertices = locs;
    surf.faces = tri;
end
    
% initialize constants and make modifications
global HL
HL = SCM_init_globs;

% [HL.Nee_a, HL.Nei_a] = deal(5000, 5000);
% lam = 1;
% HL.gamma_i = HL.gamma_i / lam;
% HL.gi = HL.gi * lam;

for k = 1:K
    
    if k <= 30 / T0 % k == 20 / T0
        source_drive = NaN;
    elseif k > 150 / T0
        source_drive = NaN;
    else
        source_drive = NaN; % 3;
    end

    fprintf(['Running simulation , ' num2str(k) ' ... \n']);
    tic; % start timer

    [samp_time,last,fine] = seizing_cortical_field(...
        source_drive, map, T0, last, ...
        locs, laplacian, avg_D, ...
        micro_idx, macro_idx, focus_idx, normal_idx, ...
        lessihb_idx);
    
    fprintf(['Run time ' num2str(toc) '\n']);

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
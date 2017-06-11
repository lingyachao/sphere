clear; close all;

%% specify run type
type = 'sphere';
note = 'bifur';
save_output = false;
visualize = true;
print_count = true;

%% load grid
if strcmp(type, 'sphere')
    load('./computed_sphere_grid/N10242_R10_wideNodes.mat'); 
    % load('./computed_sphere_grid/N42_R10.mat'); avg_D = 0.3;
    pos_hemi = locs(:,3) >= 0;
    neg_hemi = locs(:,3) < 0;
elseif strcmp(type, 'brain')
    load('./computed_brain_grid/N40962.mat');
    load('unitsphere.mat', 'coord');
else
    error('not recognized type');
end

%% set the output directory
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
if visualize
    f = figure;
    set(f, 'Position', [200 300 1300 400]);
    if strcmp(type, 'brain')
        load('autism.surface.mat', 'tri');
        surf.vertices = locs;
        surf.faces = tri;
    end
end

%% initialize parameters and map
k = 0;
K = 2000;
T0 = 0.1;
map = make_map(laplacian);

%% initialize initial state
last = make_IC(N);

%% define zones
% lessihb_idx = lessihb_area;
lessihb_idx = locs(:,3) < -6;
% lessihb_idx = coord(1,:)' > 0.5;
% lessihb_idx = true(N, 1);

zones.focus_zone = map == 1;
zones.lessihb_zone = lessihb_idx & map ~= 1;
zones.normal_zone = ~lessihb_idx;
    
%% initialize constants and make modifications
global HL
HL = SCM_init_globs(N);

% no potassium
% HL.kR = 0;

HL.k_decay = HL.k_decay * ones(N, 1);
% HL.k_decay(zones.normal_zone) = 100;

% HL.Vi_rest = -74;

last.dVe(zones.normal_zone) = -1;
% last.dVi(zones.normal_zone) = 0;
last.dVe(lessihb_idx) = -1;
% last.dVi(lessihb_idx) = 0;

last.D22(:) = 3;
% last.D11 = last.D22 / 100;

gauss_width = 0.5;
[lat,long] = GridSphere(N);
arc_dist = 10 * (lat+90) * (pi/180);
last.K = 10 * exp(-arc_dist.^2 / gauss_width.^2);

%% bifurcation stuff

% HL.ge(zones.lessihb_zone) = 0.8 * HL.ge(zones.lessihb_zone);
% HL.ge(zones.focus_zone) = 0.8 * HL.ge(zones.focus_zone);
% HL.phi_ee_sc(zones.focus_zone) = 20 * HL.phi_ee_sc(zones.focus_zone);

% gauss_width = 1;
% [lat,long] = GridSphere(N);
% arc_dist = 10 * (lat+90) * (pi/180);
% HL.phi_ee_sc = HL.phi_ee_sc(1) * (1 + 20 * exp(-arc_dist.^2 / gauss_width.^2));
% HL.ge = HL.ge(1) * (1 - 0.3 * exp(-arc_dist.^2 / gauss_width.^2));

% increase inhibitory strength in all locations other than a patch
% HL.Vi_rest(zones.normal_zone) = HL.Vi_rest(zones.normal_zone) + 3;

%% other means of introducing excitability
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

%% for plotting trace of one node
Ve_samp = [];
N_samp = 0;

%% run simulation
for k = 1:K
        
    if k == 1
        % last.K(1:7) = 10;
        % HL.Vi_rest(zones.focus_zone) = HL.Vi_rest(zones.focus_zone) + 10;
    end

    if k < 0
        source_drive = 3;
    else
        source_drive = NaN;
    end
        
    if print_count
        fprintf(['Running simulation , ' num2str(k) ' ... ']);
    end
    tic; % start timer

    [samp_time,last,fine] = seizing_cortical_field(...
        source_drive, map, T0, last, ...
        locs, laplacian, avg_D, ...
        zones, micro_idx, macro_idx, focus_idx, normal_idx);
    
    if visualize
        if strcmp(type, 'sphere')
            subplot(1, 3, 1);
            scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.K(neg_hemi), 'filled');
            % scatter(locs(pos_hemi,1), locs(pos_hemi,2), 15, last.Ve(pos_hemi), 'filled');
            % caxis([-65,-50])
            % caxis([0,30]);
            colorbar;
            title('K');

            subplot(1, 3, 2);    
            scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Ve(neg_hemi), 'filled');
            caxis([-65,-50]);
            title('Ve');
            
            subplot(1, 3, 3);    
            scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
            caxis([0,30]);
            title('Qe');
            
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
    
    % for plotting focal voltage traces
    if k == 1
        N_samp = length(fine.Ve_focus(:,1));
        Ve_samp = NaN(K * N_samp,1);
    end
    Ve_samp((k-1)*N_samp+1 : k*N_samp) = fine.Ve_focus(:,1);
    
    if print_count
        fprintf(['RT ' num2str(toc) '\n']);
        fprintf(['K normal ' num2str(mean(last.K(zones.normal_zone))) ' K abnormal ' num2str(mean(last.K(lessihb_idx))) '\n']);
        fprintf(['D2 ' num2str(mean(last.D22(lessihb_idx))) ' dVe ' num2str(mean(last.dVe(lessihb_idx))) '\n']);
        fprintf(['Ve focus ' num2str(last.Ve(1)) '\n']);
    end
end

figure;
plot(Ve_samp);

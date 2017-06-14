% load node positions on BRAIN
load('./computed_brain_grid/N40962.mat');

% load node positions on SPHERE
load('unitsphere.mat', 'coord');
locs2 = 10 * coord';
pos_hemi = locs2(:,3) >= 0;
neg_hemi = locs2(:,3) < 0;

% load faces data
surf.vertices = locs;
surf.faces = tri;
surf2.vertices = locs2;
surf2.faces = tri;
surf3 = surf2;
surf3.vertices(:,1) = -surf3.vertices(:,1);
surf3.vertices(:,3) = -surf3.vertices(:,3);

% load meta file
load(META_FILE);

% load vertices that are closest to electrodes
loc_grid_center = [64.41, -7.282, 21.48]; % center of the ECoG grid (mm)
dist_grid = 12; % distance between electrodes (mm)
[focus_idx, macro_pos, macro_idx, macro_2d, ...
            micro_pos, micro_idx, micro_2d] = ...
    generate_electrode_grid(loc_grid_center, dist_grid, RAW_DIR);
save(ELEC_FILE, 'focus_idx', 'macro_pos', 'macro_idx', 'macro_2d', 'micro_pos', 'micro_idx', 'micro_2d');

% create a filter for subsetting electrodes
[~,macro_filter] = ismember(macro_idx(:,1), lessihb_idx);
[~,micro_filter] = ismember(micro_idx, lessihb_idx);
[~,focus_filter] = ismember(focus_idx, lessihb_idx);

% allocate storing space
Qe_rand = NaN(K*T, 4);
Ve_rand = NaN(K*T, 4);
Qe_avg = NaN(K*T, 3);
Ve_avg = NaN(K*T, 3);
Qe_macro = NaN(K*T, size(macro_idx, 1));
Ve_macro = NaN(K*T, size(macro_idx, 1));
Qe_micro = NaN(K*T, size(micro_idx, 1));
Ve_micro = NaN(K*T, size(micro_idx, 1));

% start movie
vidObj = VideoWriter(VIDEO_FILE, 'MPEG-4');
vidObj.FrameRate = 23;
open(vidObj);

% set plotting window
f = figure;
set(f, 'Position', [200 300 900 400]);

for k = 1:K
    clf(f);
    
    fprintf(['Read in ' num2str(k) '\n']);
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(k) '.mat']);

    % subset macro to keep only the ones close to electrodes
    fine.Qe_focus = fine.Qe_lessihb(:,focus_filter);
    fine.Ve_focus = fine.Ve_lessihb(:,focus_filter);
    fine.Qe_macro = fine.Qe_lessihb(:,macro_filter);
    fine.Ve_macro = fine.Ve_lessihb(:,macro_filter);
    fine.Qe_micro = fine.Qe_lessihb(:,micro_filter);
    fine.Ve_micro = fine.Ve_lessihb(:,micro_filter);
    
    % save 
    Qe_rand(1+(k-1)*T : k*T,1) = fine.Qe_focus(:,1);
    Qe_rand(1+(k-1)*T : k*T,2) = fine.Qe_macro(:,1);
    Qe_rand(1+(k-1)*T : k*T,3) = fine.Qe_macro(:,25);
    Qe_rand(1+(k-1)*T : k*T,4) = fine.Qe_normal(:,1);
    Ve_rand(1+(k-1)*T : k*T,1) = fine.Ve_focus(:,1);
    Ve_rand(1+(k-1)*T : k*T,2) = fine.Ve_macro(:,1);
    Ve_rand(1+(k-1)*T : k*T,3) = fine.Ve_macro(:,25);
    Ve_rand(1+(k-1)*T : k*T,4) = fine.Ve_normal(:,1);

    Qe_avg(1+(k-1)*T : k*T,1) = fine.Qe_focus_avg;
    Qe_avg(1+(k-1)*T : k*T,2) = fine.Qe_lessihb_avg;
    Qe_avg(1+(k-1)*T : k*T,3) = fine.Qe_normal_avg;
    Ve_avg(1+(k-1)*T : k*T,1) = fine.Ve_focus_avg;
    Ve_avg(1+(k-1)*T : k*T,2) = fine.Ve_lessihb_avg;
    Ve_avg(1+(k-1)*T : k*T,3) = fine.Ve_normal_avg;

    Qe_macro(1+(k-1)*T : k*T,:) = fine.Qe_macro;
    Ve_macro(1+(k-1)*T : k*T,:) = fine.Ve_macro;
    Qe_micro(1+(k-1)*T : k*T,:) = fine.Qe_micro;
    Ve_micro(1+(k-1)*T : k*T,:) = fine.Ve_micro;
    
    % plot frame for video and write frame

    subplot(2, 3, [1,2,4,5]);
    figure_wire(surf, last.Qe, false);
    view(90, 0);
    hold on;
    % scatter3(locs(micro_idx,1), locs(micro_idx,2), locs(micro_idx,3), 15, 'g', 'filled');
    scatter3(macro_pos(:,1), macro_pos(:,2), macro_pos(:,3), 20, ...
        'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'black');
    
    subplot(2, 3, 3);
    figure_wire(surf2, last.Qe, false);
    
    subplot(2, 3, 6);
    figure_wire(surf3, last.Qe, false);
    
    drawnow;
    im = getframe(f);
    writeVideo(vidObj,im);
end

save(SAMPLE_DATA_FILE, 'Qe_rand', 'Ve_rand', 'Qe_avg', 'Ve_avg', ...
    'Qe_macro', 'Ve_macro', 'Qe_micro', 'Ve_micro');

close(vidObj);
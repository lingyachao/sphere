% load node positions on BRAIN
load('./computed_brain_grid/N40962.mat');

% load node positions on SPHERE
load('unitsphere.mat', 'coord');
locs2 = 10 * coord';
pos_hemi = locs2(:,3) >= 0;
neg_hemi = locs2(:,3) < 0;

% load vertices that are closest to electrodes
% and create a filter for subsetting macro_indices
load('./computed_brain_grid/electrode_idx.mat');
electrode_filter = ismember(macro_idx, e_ver(:,1));

% allocate storing space
Qe_rand = NaN(K*T, 3);
Ve_rand = NaN(K*T, 3);
Qe_avg = NaN(K*T, 3);
Ve_avg = NaN(K*T, 3);
Qe_macro = NaN(K*T, size(e_ver, 1));
Ve_macro = NaN(K*T, size(e_ver, 1));
Qe_micro = NaN(K*T, length(micro_idx));
Ve_micro = NaN(K*T, length(micro_idx));

% start movie
vidObj = VideoWriter(VIDEO_FILE, 'MPEG-4');
vidObj.FrameRate = 23;
open(vidObj);

% set plotting window
f = figure;
set(f, 'Position', [200 300 900 400]);
load('autism.surface.mat', 'tri');
surf.vertices = locs;
surf.faces = tri;
surf2.vertices = locs2;
surf2.faces = tri;
surf3 = surf2;
surf3.vertices(:,1) = -surf3.vertices(:,1);
surf3.vertices(:,3) = -surf3.vertices(:,3);

for k = 1:K
    clf;
    
    fprintf(['Read in ' num2str(k) '\n']);
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(k) '.mat']);

    % subset macro to keep only the ones close to electrodes
    fine.Qe_elec = fine.Qe_macro(:,electrode_filter);
    fine.Ve_macro = fine.Ve_macro(:,electrode_filter);
    
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
    scatter3(e_pos(:,1), e_pos(:,2), e_pos(:,3), 20, ...
        'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'black');
    
    subplot(2, 3, 3);
    figure_wire(surf2, last.Qe, false);
    
    subplot(2, 3, 6);
    figure_wire(surf3, last.Qe, false);
    
    drawnow;
    f = getframe;
    writeVideo(vidObj,f);
end

save(SAMPLE_DATA_FILE, 'Qe_rand', 'Ve_rand', 'Qe_avg', 'Ve_avg', ...
    'Qe_macro', 'Ve_macro', 'Qe_micro', 'Ve_micro');

close(vidObj);
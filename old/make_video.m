% specifiy directory where simulation results were saved
DATA_DIR = './data/sphere_N10242_R10_04251642_normal_1.2_stimulus_3/';

% load node positions
load('N10242_R10_wideNodes.mat');
% load('unitsphere.mat', 'locs');
% locs = 10 * coord';
pos_hemi = locs(:,3) >= 0;
neg_hemi = locs(:,3) < 0;

% start movie
vidObj = VideoWriter([DATA_DIR 'movie_sparse_sphere.mp4'], 'MPEG-4');
vidObj.FrameRate=23;
open(vidObj);

% set plotting window
f = figure;
set(f, 'Position', [200 300 900 400]);

for i = 1:2000
    fprintf(['Read in ' num2str(i) '\n']);    
    load([DATA_DIR 'seizing_cortical_field_k_'  num2str(i) '.mat'], 'last');

    % subplot(1, 2, 1);
    scatter(locs(pos_hemi,1), locs(pos_hemi,2), 15, last.Qe(pos_hemi), 'filled');
    caxis([0,30]);
    xlim([-10,33]);
    
    % subplot(1, 2, 2);
    hold on;
    scatter(23 - locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
    caxis([0,30]);
    drawnow;
        
    % scatter(23 - locs(micro_idx,1),locs(micro_idx,2), 15, 'g', 'filled');
    scatter(23 - locs(macro_idx,1), locs(macro_idx,2), 15, 'r', 'filled');
    
    f = getframe;
    writeVideo(vidObj,f);
    clf;
end

close(vidObj);
% specifiy directory where simulation results were saved
DATA_DIR = './data/brain_N40962_03071708/';

% choose the map type, and specify a starting index to load
map_type = 'fixed_point_source';

% load node positions on BRAIN
load('./computed_brain_grid/N40962.mat');

% load node positions on SPHERE
load('unitsphere.mat', 'coord');
locs2 = 10 * coord';
pos_hemi = locs2(:,3) >= 0;
neg_hemi = locs2(:,3) < 0;

% start movie
vidObj = VideoWriter([DATA_DIR 'movie_sparse_combine.avi']);
vidObj.FrameRate = 1;
open(vidObj);

% set plotting window
f = figure;
set(f, 'Position', [200 0 1200 800]);
load('autism.surface.mat', 'tri');
surf.vertices = locs;
surf.faces = tri;
surf2.vertices = locs2;
surf2.faces = tri;
surf3 = surf2;
surf3.vertices(:,1) = -surf3.vertices(:,1);
surf3.vertices(:,3) = -surf3.vertices(:,3);

for i = 1:200
    fprintf(['Read in ' num2str(i) '\n']);    
    load([DATA_DIR 'seizing_cortical_field_k_'  num2str(i) '.mat'], 'Qe');

    clf(f);
    
    subplot(2, 3, [1,2,4,5]);
    figure_wire(surf, locs(:,3)); 
    
    subplot(2, 3, 3);
    figure_wire(surf2, Qe);
    
    subplot(2, 3, 6);
    figure_wire(surf3, Qe);
    drawnow;
    
    im = getframe(f);
    writeVideo(vidObj,im);
end

close(vidObj);
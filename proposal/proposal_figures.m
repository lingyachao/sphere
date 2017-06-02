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

  
load([DATA_DIR 'seizing_cortical_field_k_2.mat'], 'Qe');

subplot(1,2,1);
figure_wire(surf, Qe, true);

subplot(1,2,2);
figure_wire(surf, Qe, false);

%% 


% load node positions on BRAIN
load('./computed_sphere_grid/N10242_R10.mat');

% set plotting window
f = figure;
set(f, 'Position', [200 0 1200 800]);
load('autism.surface.mat', 'tri');
surf.vertices = locs;
surf.faces = tri;

figure_wire(surf, Qe, true);

%%


[V,F] = icosphere(3);
surf.vertices = V;
surf.faces = F;
figure_wire(surf, NaN, true);


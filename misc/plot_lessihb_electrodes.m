load('N40962.mat'); avg_D = 0.3777;
load('unitsphere.mat');

surf.vertices = locs;
surf.faces = tri;
surf_sphere.vertices = 10 * coord';
surf_sphere.faces = tri;

col = zeros(N, 1);
col(lessihb_area_brain) = 30;

figure_wire(surf, col, false);
view(90, 0);
hold on;

scatter3(macro_pos(:,1), macro_pos(:,2), macro_pos(:,3), 20, ...
    'filled', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'black');
scatter3(micro_pos(:,1), micro_pos(:,2), micro_pos(:,3), 20, ...
    'filled', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'black');
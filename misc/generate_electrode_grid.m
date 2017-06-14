% load vertice locations
load('./computed_brain_grid/N40962.mat', 'N', 'locs');
load('autism.surface.mat', 'tri');
surf.vertices = locs;
surf.faces = tri;

% figure_wire(surf, NaN, true);

% get boundary of brain
k = boundary(locs, 0.3);
pt = trisurf(k, locs(:,1), locs(:,2), locs(:,3),'Facecolor','w');

% compute normal at all vertices
TR = triangulation(pt.Faces, pt.Vertices);
VN = vertexNormal(TR);

% choose a center electrode location
selected = [64.41, -7.282, 21.48];

% find a patch or nearby vertices and average the norm
to_avg = find(abs(locs(:,1) - selected(1)) < 0.01 & abs(locs(:,2) - selected(2)) < 0.01 ...
    & abs(locs(:,3) - selected(3)) < 0.01);

all_pairs = [pt.Faces(:,[1,2]); ...
             pt.Faces(:,[1,3]); ...
             pt.Faces(:,[2,1]); ...
             pt.Faces(:,[2,3]); ...
             pt.Faces(:,[3,1]); ...
             pt.Faces(:,[3,2])];

for itr = 1:3
    for c = to_avg'
        to_avg = union(to_avg, all_pairs(all_pairs(:,1) == c,2));
    end
end

n = mean(VN(to_avg,:));
hold on;
quiver3(selected(1), selected(2), selected(3), n(1), n(2), n(3), 15, 'Color', 'r');

%% get electrode locations

clf;
load('./data/brain_N40962_06131715_full_data/raw/seizing_cortical_field_k_50.mat');
figure_wire(surf, last.Qe, false);
hold on;

perp = null(n)';
dist_grid = 12;

quiver3(selected(1), selected(2), selected(3), n(1), n(2), n(3), 15, 'Color', 'red');

center = selected + 0.5*n;
e_pos = ones(25, 1) * center + ...
    dist_grid * repelem(-2:2, 5)' * perp(1,:) + ...
    dist_grid * repmat((-2:2)', 5, 1) * perp(2,:);

scatter3(e_pos(:,1), e_pos(:,2), e_pos(:,3), 40, ...
    'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'black');

%% find nearest vertices for each electrode

e_ver = NaN(25, 5);

for k = 1:25
    dist = locs - ones(N, 1)*e_pos(k,:);
    dist = sum(abs(dist).^2, 2).^(1/2);
    [~,I] = sort(dist);
    e_ver(k,:) = I(1:5);
end

cols = {'r', 'm', 'g'};
for k = 1:3
    scatter3(locs(e_ver(:,k),1), locs(e_ver(:,k),2), locs(e_ver(:,k),3), 40, ...
        'filled', 'MarkerFaceColor', cols{k}, 'MarkerEdgeColor', 'black');
end

save('./computed_brain_grid/electrode_idx.mat', 'e_ver', 'e_pos');











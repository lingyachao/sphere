function [focus_idx, macro_pos, macro_idx, macro_2d, ...
                     micro_pos, micro_idx, micro_2d] = generate_electrode_grid( ...
    loc_grid_center, dist_grid, RAW_DIR)

    load('./computed_brain_grid/N40962.mat');

    % get boundary of brain
    k = boundary(locs, 0.3);
    pt = trisurf(k, locs(:,1), locs(:,2), locs(:,3),'Facecolor','w');

    % compute normal at all vertices
    TR = triangulation(pt.Faces, pt.Vertices);
    VN = vertexNormal(TR);

    % find a patch or nearby vertices and average the norm
    to_avg = find(abs(locs(:,1) - loc_grid_center(1)) < 0.01 & abs(locs(:,2) - loc_grid_center(2)) < 0.01 ...
        & abs(locs(:,3) - loc_grid_center(3)) < 0.01);

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

    %% get electrode locations

    clf;
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(50) '.mat']);
    surf.vertices = locs;
    surf.faces = tri;
    figure_wire(surf, last.Qe, false);
    hold on;
    quiver3(loc_grid_center(1), loc_grid_center(2), loc_grid_center(3), n(1), n(2), n(3), 15, 'Color', 'red');

    center = loc_grid_center + 0.5*n;
    perp = null(n)';
    macro_pos = ones(25, 1) * center + ...
        dist_grid * repelem(-2:2, 5)' * perp(1,:) + ...
        dist_grid * repmat((-2:2)', 5, 1) * perp(2,:);
    macro_2d = [repelem(-2:2, 5)' repmat((-2:2)', 5, 1)];

    scatter3(macro_pos(:,1), macro_pos(:,2), macro_pos(:,3), 40, ...
        'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'black');

    %% find nearest vertices for each electrode

    macro_idx = NaN(25, 5);

    for k = 1:25
        dist = locs - ones(N, 1)*macro_pos(k,:);
        dist = sum(abs(dist).^2, 2).^(1/2);
        [~,I] = sort(dist);
        macro_idx(k,:) = I(1:5);
    end

    cols = {'r', 'm', 'g'};
    for k = 1:3
        scatter3(locs(macro_idx(:,k),1), locs(macro_idx(:,k),2), locs(macro_idx(:,k),3), 40, ...
            'filled', 'MarkerFaceColor', cols{k}, 'MarkerEdgeColor', 'black');
    end

    %% save focus, macro and micro data
    focus_idx = find(make_map(laplacian));
    micro_pos = [];
    micro_idx = [];
    micro_2d = [];
    
end











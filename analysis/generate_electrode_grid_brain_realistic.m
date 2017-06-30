function [focus_idx, macro_pos, macro_idx, macro_2d, ...
                     micro_pos, micro_idx, micro_2d] = generate_electrode_grid_brain_realistic( ...
    loc_grid_center, dist_grid, RAW_DIR)

    % loc_grid_center = [64.41, -7.28, 21.48];
    % dist_grid = 12;
    % RAW_DIR = './data/brain_N40962_06141205_full_data/raw/';

    load('N40962.mat');
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(50) '.mat']);

    % get boundary of brain
    k = boundary(locs, 0.3);
    pt = trisurf(k, locs(:,1), locs(:,2), locs(:,3), last.Qe);
    view(90, 0);

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
    surf.vertices = locs;
    surf.faces = tri;
    figure_wire(surf, last.Qe, false);
    view(90, 0);
    hold on;
    quiver3(loc_grid_center(1), loc_grid_center(2), loc_grid_center(3), n(1), n(2), n(3), 15, 'Color', 'red');

    center = loc_grid_center + 0*n;
    perp = null(n)';
    macro_pos = ones(25, 1) * center + ...
        dist_grid * repelem(-2:2, 5)' * perp(1,:) + ...
        dist_grid * repmat((-2:2)', 5, 1) * perp(2,:);
    macro_2d = [repelem(-2:2, 5)' repmat((-2:2)', 5, 1)];

    scatter3(macro_pos(:,1), macro_pos(:,2), macro_pos(:,3), 40, ...
        'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'black');

    %% find nearest vertices for each electrode

    macro_idx = zeros(25, N);

    VN2 = vertexNormal(triangulation(tri, locs));

    for k = 1:25
        ur = ones(N, 1) * macro_pos(k,:) - locs;
        dist = sum(ur.^2, 2).^(1/2);
        ur = ur ./ (dist * ones(1, 3));
        pot_coeff = sum(ur.*VN2, 2) ./ dist.^2;
        % scatter3(locs(I(1:10),1), locs(I(1:10),2), locs(I(1:10),3), 40, ...
        %     'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'black');

        [~,I] = sort(dist);
        macro_idx(k,I(1:25)) = pot_coeff(I(1:25));
    end

    % cols = {'r', 'm', 'g'};
    % for k = 1:3
    %     scatter3(locs(macro_proj(:,k),1), locs(macro_proj(:,k),2), locs(macro_proj(:,k),3), 15, ...
    %         'filled', 'MarkerFaceColor', cols{k}, 'MarkerEdgeColor', 'black');
    % end

    %% save focus, macro and micro data
    focus_idx = find(make_map(laplacian));
    micro_pos = [];
    micro_idx = [];
    micro_2d = [];
    
end
function [focus_idx, macro_pos, macro_transform, macro_2d, ...
                     micro_pos, micro_transform, micro_2d] = generate_electrode_grid_brain( ...
    loc_grid_center, dist_grid, RAW_DIR, flag_dipole)

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

    center = loc_grid_center + 2*n;
    perp = null(n)';
    macro_pos = ones(25, 1) * center + ...
        dist_grid * repelem(-2:2, 5)' * perp(1,:) + ...
        dist_grid * repmat((-2:2)', 5, 1) * perp(2,:);
    macro_2d = [repelem(-2:2, 5)' repmat((-2:2)', 5, 1)];

    scatter3(macro_pos(:,1), macro_pos(:,2), macro_pos(:,3), 40, ...
        'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'black');

    %% find nearest vertices for each electrode

    macro_transform = zeros(25, N);
    
    
    for k = 1:25
        ur = ones(N, 1) * macro_pos(k,:) - locs;
        dist = sum(ur.^2, 2).^(1/2);
        ur = ur ./ (dist * ones(1, 3));
        [~,I] = sort(dist);
        
        if flag_dipole
            VN2 = vertexNormal(triangulation(tri, locs));
            pot_coeff = sum(ur.*VN2, 2) ./ dist.^2;
            macro_transform(k,I(1:25)) = pot_coeff(I(1:25));
        else
            top_N = 7;
            macro_transform(k,I(1:top_N)) = 1/top_N;

            cols = {'r', 'm', 'g'};
            for close = 1:3
                scatter3(locs(I(close),1), locs(I(close),2), locs(I(close),3), 15, ...
                    'filled', 'MarkerFaceColor', cols{close}, 'MarkerEdgeColor', 'black');
            end
        end
    end
    
    %% save focus, macro and micro data
    focus_idx = find(make_map(laplacian));
    micro_pos = zeros(0, 3);
    micro_transform = zeros(0, N);
    micro_2d = zeros(0, 2);
    
end











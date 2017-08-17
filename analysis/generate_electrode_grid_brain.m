function [focus_idx,macro_pos,macro_transform,macro_2d, ...
                    micro_pos,micro_transform,micro_2d] = generate_electrode_grid_brain( ...
    loc_grid_center, dist_grid_macro, dist_grid_micro, RAW_DIR, flag_dipole, closest_N)

    load('N40962.mat');
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(2000) '.mat']);

    % get boundary of brain
    figure;
    bounds = boundary(locs, 0.3);
    pt = trisurf(bounds, locs(:,1), locs(:,2), locs(:,3), last.Qe);
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

    ctr_props.center = loc_grid_center + 1*n;
    ctr_props.perp = null(n)';
    ctr_props.grid_basis = [repelem(-2:2, 5)' repmat((-2:2)', 5, 1)];
    
    [macro_pos,macro_transform,macro_2d] = gain_matrix(N, locs, dist_grid_macro, ctr_props, flag_dipole, closest_N);
    [micro_pos,micro_transform,micro_2d] = gain_matrix(N, locs, dist_grid_micro, ctr_props, flag_dipole, closest_N);
    focus_idx = find(make_map(laplacian));
    
    % micro_pos = zeros(0, 3);
    % micro_transform = zeros(0, N);
    % micro_2d = zeros(0, 2);
    
end

function [pos,transform,pos_2d] = gain_matrix(N, locs, dist_grid, ctr_props, flag_dipole, closest_N)

    pos = ones(25, 1) * ctr_props.center + ...
        dist_grid * repelem(-2:2, 5)' * ctr_props.perp(1,:) + ...
        dist_grid * repmat((-2:2)', 5, 1) * ctr_props.perp(2,:);
    pos_2d =  dist_grid * ctr_props.grid_basis;

    scatter3(pos(:,1), pos(:,2), pos(:,3), 40, ...
        'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'black');

    transform = zeros(25, N);

    for k = 1:25
        ur = ones(N, 1) * pos(k,:) - locs;
        dist = sum(ur.^2, 2).^(1/2);
        % ur = ur ./ (dist * ones(1, 3));
        [~,I] = sort(dist);
        
        if flag_dipole
            % VN2 = vertexNormal(triangulation(tri, locs));
            % pot_coeff = sum(ur.*VN2, 2) ./ dist.^2;
            pot_coeff = ones(N, 1) ./ dist;
            top_coeff = pot_coeff(I(1:closest_N));
            transform(k,I(1:closest_N)) = top_coeff / sum(top_coeff);
        else
            transform(k,I(1:closest_N)) = 1/closest_N;
        end
        
%         cols = {'r', 'm', 'g'};
%         for close = 1:3
%             scatter3(locs(I(close),1), locs(I(close),2), locs(I(close),3), 15, ...
%                 'filled', 'MarkerFaceColor', cols{close}, 'MarkerEdgeColor', 'black');
%         end        
    end
end











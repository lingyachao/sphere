function find_electrode_grid_center(RAW_DIR, dist_grid)

    load('N40962.mat', 'N', 'locs', 'tri');
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(50) '.mat'], 'last');

    % get boundary of brain
    bounds = boundary(locs, 0.3);
    pt = trisurf(bounds, locs(:,1), locs(:,2), locs(:,3), last.Qe);
    
    fig = figure;
    set(fig, 'Position', [200 300 800 400]);
    surf.vertices = locs;
    surf.faces = tri;
    figure_wire(surf, last.Qe, false);
    view(90, 0);
    hold on;

    dcm_obj = datacursormode(fig);
    set(dcm_obj, 'UpdateFcn', @plot_grid);
    set(dcm_obj, 'enable', 'on');
    
    function txt = plot_grid(~, event_obj)
        txt = '';
        
        loc_grid_center = get(event_obj,'Position');
        
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

        % get electrode locations

        % close(fig);
        % s = subplot(1,2,2); cla(s);
        % fig = figure;
        % set(fig, 'Position', [200 300 800 400]);
        % figure_wire(surf, last.Qe, false);
        % view(90, 0);
        quiver3(loc_grid_center(1), loc_grid_center(2), loc_grid_center(3), n(1), n(2), n(3), 15, 'Color', 'red');

        center = loc_grid_center + 1*n;
        perp = null(n)';
        macro_pos = ones(25, 1) * center + ...
            dist_grid * repelem(-2:2, 5)' * perp(1,:) + ...
            dist_grid * repmat((-2:2)', 5, 1) * perp(2,:);

        scatter3(macro_pos(:,1), macro_pos(:,2), macro_pos(:,3), 40, ...
            'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'black');
        
        % print info
        fprintf(['[' num2str(loc_grid_center(1), '%.2f') ', ' ...
                     num2str(loc_grid_center(2), '%.2f') ', '... 
                     num2str(loc_grid_center(3), '%.2f') '] \n'])
    end
end











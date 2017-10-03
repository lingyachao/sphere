function selected = interactive_find_node_boundary_brain
    clear; close all;

    load('./computed_brain_grid/N40962.mat', 'N', 'locs');
    load('./computed_sphere_grid/unitsphere.mat', 'coord');
    % [lat,~] = GridSphere(10242);

    p = coord(1,:)' > 0;
    fig = figure;
    
    idxs = 1:N;
    idxs = idxs(p);
    colors = ones(N, 1);
    colors(coord(1,:)' > 0.7) = 2;
    % colors(macro_idx) = 0;
    s = scatter(coord(2,p)', coord(3,p)', 15, colors(p), 'filled');
    hold on;
    
    selected = [];
    
    dcm_obj = datacursormode(fig);
    set(dcm_obj, 'UpdateFcn', @give_idx);
    set(dcm_obj, 'enable', 'on');
    
    function txt = give_idx(~, event_obj)
        I = get(event_obj, 'DataIndex');
        I = idxs(I);
        
        disp(['chose' num2str(I)]);

        selected = [selected, I];
        colors(I) = 0;
        selected = sort(unique(selected));
        set(s, 'CData', colors(p));
        
        disp(mat2str(selected));
        txt = {''};
    end
end
function selected = find_node_pos_ihb_boundary
    clear; close all;

    load('./computed_sphere_grid/N10242_R10.mat', 'N', 'locs');
    [lat,~] = GridSphere(10242);

    p = lat < 0;
    fig = figure;
    
    idxs = 1:N;
    idxs = idxs(p);
    colors = ones(N, 1);
    % colors(macro_idx) = 0;
    s = scatter(locs(p,1),locs(p,2), 15, colors(p), 'filled');
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
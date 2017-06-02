function selected = find_node_pos
    % clear; close all;

    load('./computed_sphere_grid/N10242_R10.mat', 'N', 'locs');
    [lat,long] = GridSphere(10242);

    p = lat < -40 & lat > -80 & long > 0;
    fig = figure;
    
    idxs = 1:N;
    idxs = idxs(p);
    colors = ones(N, 1);
    s = scatter(locs(p,1),locs(p,2), 15, colors(p), 'filled');
    hold on;
    
    t = (1:10242)';
    labels = num2str(t(p),'%d'); 
    text(locs(p,1),locs(p,2),  labels, 'horizontal', 'left', 'vertical', 'bottom', 'FontSize', 5);

    selected = [];
    
    dcm_obj = datacursormode(fig);
    set(dcm_obj, 'UpdateFcn', @give_idx);
    set(dcm_obj, 'enable', 'on');
    
    function txt = give_idx(~, event_obj)
        I = get(event_obj, 'DataIndex');
        I = idxs(I);
        
        disp(['chose' num2str(I)]);
        
        if ismember(I, selected)
            selected = selected(selected ~= I);
            colors(I) = 1;
        else
            selected = [selected, I];
            colors(I) = 0;
        end
        % selected = sort(selected);
        set(s, 'CData', colors(p));
        
        disp(mat2str(selected));
        txt = {''};
    end
end

% hold on;
% micro_idx = [5121,5311,5312,5120,5122,4931,4932];
% micro_idx = [247, ...
%              303, 287, 233, 197, 230, 284, ...
%              361, 343, 335, 275, 217, 181, 153, 178, 214, 272, 332, 340];
% scatter3(locs(micro_idx,1),locs(micro_idx,2),locs(micro_idx,3), 15, 'b', 'filled');
% 
% % macro_idx = [5121,5803,5804,5117,5125,4439,4440,5493,5494,5119,5123,4749,4750];
% macro_idx = [247, ...
%              361, 335, 217, 153, 214, 332, ...
%              495, 465, 445, 315, 213, 131, 79, 128, 210, 312, 442, 462];
% scatter3(locs(macro_idx,1),locs(macro_idx,2),locs(macro_idx,3), 15, 'r', 'filled');
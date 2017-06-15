function plot_brain_instance(surf, surf_sphere, last, macro_pos, micro_pos)
    
    subplot(2, 3, [1,2,4,5]);
    figure_wire(surf, last.Qe, false);
    view(90, 0);
    hold on;
    
    if ~isnan(macro_pos)
        scatter3(macro_pos(:,1), macro_pos(:,2), macro_pos(:,3), 20, ...
            'filled', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'black');
    end
    if ~isnan(micro_pos)
        scatter3(micro_pos(:,1), micro_pos(:,2), micro_pos(:,3), 20, ...
            'filled', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'black');
    end
    
    subplot(2, 3, 3);
    figure_wire(surf_sphere, last.Qe, false);
    view(90, 0);
    
    subplot(2, 3, 6);
    figure_wire(surf_sphere, last.Qe, false);
    view(270, 0);
    
end


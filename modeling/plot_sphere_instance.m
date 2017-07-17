function plot_sphere_instance(locs, last, macro_pos, micro_pos)
    
    neg_hemi = locs(:,3) < 0;

    subplot(1, 2, 1);
    scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.K(neg_hemi), 'filled');
    caxis([0,1]); axis off;

    subplot(1, 2, 2);    
    scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
    caxis([0,30]); axis off;
    
    hold on;
    if ~isnan(macro_pos)
        scatter(macro_pos(:,1), macro_pos(:,2), 15, ...
            'filled', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'black');
    end
    if ~isnan(micro_pos)
        scatter(micro_pos(:,1), micro_pos(:,2), 15, ...
            'filled', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'black');
    end
    
end


function plot_sphere_instance(locs, last, macro_idx, micro_idx)
    
    neg_hemi = locs(:,3) < 0;

    subplot(1, 3, 1);
    scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.K(neg_hemi), 'filled');
    caxis([0,20]); axis off;

    subplot(1, 3, 2);    
    scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
    caxis([0,30]); axis off;
    
    subplot(1, 3, 3);    
    scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qi_fs(neg_hemi), 'filled');
    caxis([0,120]); axis off;
    
    hold on;
    if ~isnan(macro_idx)
        scatter(locs(macro_idx,1), locs(macro_idx,2), 15, ...
            'filled', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'black');
    end
    if ~isnan(micro_idx)
        scatter(locs(micro_idx,1), locs(micro_idx,2), 15, ...
            'filled', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'black');
    end
    
end


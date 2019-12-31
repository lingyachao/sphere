% specify time points to plot
ks = int32((10:10:total_time) / T0);
ncols = ceil(length(ks) / 4);

for i = 1:length(ks)                                         
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(ks(i)) '.mat']);
    
    figure(fg_course_neg); axes('parent', tab_course_neg);
    subplot(4, ncols, i);
    
    if strcmp(type, 'sphere')
        scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
        caxis([0,30]); axis off; hold on;
        scatter(macro_pos(:,1), macro_pos(:,2), 10, 'filled', 'MarkerFaceColor', 'm');
    else
        figure_wire(surf, last.Qe, false);
        view(90, 0);
    end
    title(['t = ' num2str(ks(i) * T0)]);
    
    figure(fg_course_pos); axes('parent', tab_course_pos);
    subplot(4, ncols, i);
    
    if strcmp(type, 'sphere')
        scatter(locs(pos_hemi,1), locs(pos_hemi,2), 15, last.Qe(pos_hemi), 'filled');
        caxis([0,30]); axis off;
        title(['t = ' num2str(ks(i) * T0)]);
    end
    
    
end
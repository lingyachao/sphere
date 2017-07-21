fg_course = figure;

% specify time points to plot
ks = int32(K * (1/P:1/P:1));
ncols = ceil(length(ks) / 4);

for i = 1:length(ks)                                      
    % fprintf(['Read in ' num2str(ks(i)) '\n']);    
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(ks(i)) '.mat']);
    subplot(4, ncols, i);

    if strcmp(type, 'sphere')
        scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
        caxis([0,30]); axis off;
        hold on;
        scatter(macro_pos(:,1), macro_pos(:,2), 10, 'filled', 'MarkerFaceColor', 'm');
    else
        figure_wire(surf_sphere, last.Qe, false);
        view(90, 0);
    end
    
    title(['t = ' num2str(ks(i) * T0)]);
end

% save figure
saveas(fg_course, COURSE_FIG);
figure;

% specify time points to plot
ks = K/20:K/20:K;
ncols = ceil(length(ks) / 2);

% load node positions
load('./computed_sphere_grid/N10242_R10.mat');
pos_hemi = locs(:,3) >= 0;
neg_hemi = locs(:,3) < 0;

for i = 1:length(ks)                                       
    fprintf(['Read in ' num2str(ks(i)) '\n']);    
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(ks(i)) '.mat']);

    if i <= ncols
        subplot(5, ncols, i);
    else
        subplot(5, ncols, i + 2*ncols)
    end
    scatter(locs(pos_hemi,1), locs(pos_hemi,2), 15, last.Qe(pos_hemi), 'filled');
    title(['t = ' num2str(ks(i) * T0)])
    caxis([0,30]); axis off;
    
    if i <= ncols
        subplot(5, ncols, i + ncols); 
    else
        subplot(5, ncols, i + 3*ncols);
    end
    scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
    caxis([0,30]); axis off;    
    drawnow;
end
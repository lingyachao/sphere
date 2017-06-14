load('./data/sphere_N40962_R10_02241642/seizing_cortical_field_k_11.mat');
load('./sphere_grid/N40962_R10.mat');

pos_hemi = locs(:,3) >= 0;
neg_hemi = locs(:,3) < 0;

f = figure;
set(f, 'Position', [200 300 900 400]);

subplot(1, 2, 1);
scatter(locs(pos_hemi,1), locs(pos_hemi,2), 15, Qe(pos_hemi), 'filled');
caxis([0,30]);
subplot(1, 2, 2);    
scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, Qe(neg_hemi), 'filled');
caxis([0,30]);
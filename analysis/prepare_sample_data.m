% create a filter for subsetting focal nodes
[~,focus_filter] = ismember(focus_idx, lessihb_idx);

% allocate storing space
Qe_rand  = NaN(K*T, 4);
Ve_rand  = NaN(K*T, 4);
Qe_avg   = NaN(K*T, 3);
Ve_avg   = NaN(K*T, 3);

Qe_macro = NaN(K*T, size(macro_transform, 1));
Ve_macro = NaN(K*T, size(macro_transform, 1));
Qe_micro = NaN(K*T, size(micro_transform, 1));
Ve_micro = NaN(K*T, size(micro_transform, 1));

% single node
node_id = 500;
[Qe_1, Qi_1, Ve_1, Vi_1, D22_1, dVe_1, dVi_1, K_1] = deal(NaN(K, 1));

% start movie
vidObj = VideoWriter(VIDEO_FILE, 'MPEG-4');
vidObj.FrameRate = 23;
open(vidObj);

% set plotting window
f = figure;
set(f, 'Position', [200 300 900 400]);

for k = 1:K
    
    fprintf(['Read in ' num2str(k) '\n']);
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(k) '.mat']);

    % subset macro to keep only the ones close to electrodes
    fine.Qe_focus = fine.Qe_lessihb(:,focus_filter);
    fine.Ve_focus = fine.Ve_lessihb(:,focus_filter);
    fine.Qe_macro = fine.Qe_lessihb * macro_transform(:,lessihb_idx)';
    fine.Ve_macro = fine.Ve_lessihb * macro_transform(:,lessihb_idx)';
    fine.Qe_micro = fine.Qe_lessihb * micro_transform(:,lessihb_idx)';
    fine.Ve_micro = fine.Ve_lessihb * micro_transform(:,lessihb_idx)';
    
    % fill in data 
    Qe_rand(1+(k-1)*T : k*T,1) = fine.Qe_focus(:,1);
    Qe_rand(1+(k-1)*T : k*T,2) = fine.Qe_macro(:,1);
    Qe_rand(1+(k-1)*T : k*T,3) = fine.Qe_macro(:,8);
    % Qe_rand(1+(k-1)*T : k*T,4) = fine.Qe_normal(:,1);
    Ve_rand(1+(k-1)*T : k*T,1) = fine.Ve_focus(:,1);
    Ve_rand(1+(k-1)*T : k*T,2) = fine.Ve_macro(:,1);
    Ve_rand(1+(k-1)*T : k*T,3) = fine.Ve_macro(:,8);
    % Ve_rand(1+(k-1)*T : k*T,4) = fine.Ve_normal(:,1);

    Qe_avg(1+(k-1)*T : k*T,1) = fine.Qe_focus_avg;
    Qe_avg(1+(k-1)*T : k*T,2) = fine.Qe_lessihb_avg;
    Qe_avg(1+(k-1)*T : k*T,3) = fine.Qe_normal_avg;
    Ve_avg(1+(k-1)*T : k*T,1) = fine.Ve_focus_avg;
    Ve_avg(1+(k-1)*T : k*T,2) = fine.Ve_lessihb_avg;
    Ve_avg(1+(k-1)*T : k*T,3) = fine.Ve_normal_avg;

    Qe_macro(1+(k-1)*T : k*T,:) = fine.Qe_macro;
    Ve_macro(1+(k-1)*T : k*T,:) = fine.Ve_macro;
    Qe_micro(1+(k-1)*T : k*T,:) = fine.Qe_micro;
    Ve_micro(1+(k-1)*T : k*T,:) = fine.Ve_micro;
    
    % single node
    Qe_1(k) = last.Qe(node_id);
    Qi_1(k) = last.Qi(node_id);
    Ve_1(k) = last.Ve(node_id);
    Vi_1(k) = last.Vi(node_id);
    D22_1(k) = last.D22(node_id);
    dVe_1(k) = last.dVe(node_id);
    dVi_1(k) = last.dVi(node_id);
    K_1(k) = last.K(node_id);
    
    % plot frame for video and write frame
    clf(f);
    if strcmp(type, 'sphere')
        plot_sphere_instance(locs, last, macro_pos, micro_pos);
    else
        plot_brain_instance(surf, surf_sphere, last, macro_pos, micro_pos);
    end
    drawnow;
    im = getframe(f);
    writeVideo(vidObj,im);
    
%     for i = 1:T
%         clf(f);
%         
%         subplot(1, 2, 1);
%         scatter(locs(pos_hemi,1), locs(pos_hemi,2), 15, fine.Qe_lessihb(i, pos_hemi)', 'filled');
%         caxis([0,30]); axis off;
% 
%         subplot(1, 2, 2);    
%         scatter(locs(neg_hemi,1), locs(neg_hemi,2), 15, fine.Qe_lessihb(i, neg_hemi)', 'filled');
%         caxis([0,30]); axis off;
%         
%         drawnow;
%         im = getframe(f);
%         writeVideo(vidObj,im);
%     end
end

single_node = table(Qe_1, Qi_1, Ve_1, Vi_1, D22_1, dVe_1, dVi_1, K_1);

save(SAMPLE_DATA_FILE, ...
    'Qe_rand', 'Ve_rand', 'Qe_avg', 'Ve_avg', ...
    'Qe_macro', 'Ve_macro', 'Qe_micro', 'Ve_micro', 'single_node');

close(vidObj);
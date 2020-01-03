% create a filter for subsetting focal nodes
[~,focus_filter] = ismember(focus_idx, fine_idx);

% allocate storing space
Qe_rand  = NaN(K*T, 3);
Ve_rand  = NaN(K*T, 3);
Qe_avg   = NaN(K*T, 3);
Ve_avg   = NaN(K*T, 3);

Qe_macro = NaN(K*T, size(macro_transform, 1));
Ve_macro = NaN(K*T, size(macro_transform, 1));
Qe_micro = NaN(K*T, size(micro_transform, 1));
Ve_micro = NaN(K*T, size(micro_transform, 1));

if exist('fine.sn_Qe', 'var')
    [sn.Qe, sn.Ve, sn.dVe, ...
        sn.Qi, sn.Vi, sn.dVi, ...
        sn.Qi_fs, sn.Vi_fs, sn.dVi_fs, ...
        sn.D22, sn.K] = deal(NaN(K*T, 1));
else
    [sn.Qe, sn.Ve, sn.dVe, ...
        sn.Qi, sn.Vi, sn.dVi, ...
        sn.Qi_fs, sn.Vi_fs, sn.dVi_fs, ...
        sn.D22, sn.K] = deal(NaN(K, 1));
end

if flag_video
    % start movie
    vidObj = VideoWriter(VIDEO_FILE, 'MPEG-4');
    vidObj.FrameRate = 23;
    open(vidObj);

    % set plotting window
    f = figure;
    set(f, 'Position', [200 300 900 400]);
end

for k = 1:K
    
    % fprintf(['Read in ' num2str(k) '\n']);
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(k) '.mat']);

    % subset macro to keep only the ones close to electrodes
    fine.Qe_focus = fine.Qe_lessihb(:,focus_filter);
    fine.Ve_focus = fine.Ve_lessihb(:,focus_filter);
    fine.Qe_macro = fine.Qe_lessihb * macro_transform(:,fine_idx)';
    fine.Ve_macro = fine.Ve_lessihb * macro_transform(:,fine_idx)';
    fine.Qe_micro = fine.Qe_lessihb * micro_transform(:,fine_idx)';
    fine.Ve_micro = fine.Ve_lessihb * micro_transform(:,fine_idx)';
    
    % fill in data 
    Qe_rand(1+(k-1)*T : k*T,1) = fine.Qe_focus(:,1);
    Qe_rand(1+(k-1)*T : k*T,2) = fine.Qe_macro(:,1);
    Qe_rand(1+(k-1)*T : k*T,3) = fine.Qe_macro(:,13);
    Ve_rand(1+(k-1)*T : k*T,1) = fine.Ve_focus(:,1);
    Ve_rand(1+(k-1)*T : k*T,2) = fine.Ve_macro(:,1);
    Ve_rand(1+(k-1)*T : k*T,3) = fine.Ve_macro(:,13);

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
    if exist('fine.sn_Qe', 'var')
        sn.Qe(1+(k-1)*T : k*T,:) = fine.sn_Qe;
        sn.Ve(1+(k-1)*T : k*T,:) = fine.sn_Ve;
        sn.dVe(1+(k-1)*T : k*T,:) = fine.sn_dVe;
        sn.Qi(1+(k-1)*T : k*T,:) = fine.sn_Qi;
        sn.Vi(1+(k-1)*T : k*T,:) = fine.sn_Vi;
        sn.dVi(1+(k-1)*T : k*T,:) = fine.sn_dVi;
        sn.Qi_fs(1+(k-1)*T : k*T,:) = fine.sn_Qi_fs;
        sn.Vi_fs(1+(k-1)*T : k*T,:) = fine.sn_Vi_fs;
        sn.dVi_fs(1+(k-1)*T : k*T,:) = fine.sn_dVi_fs;
        sn.D22(1+(k-1)*T : k*T,:) = fine.sn_D22;
        sn.K(1+(k-1)*T : k*T,:) = fine.sn_K;
    else
        if strcmp(type, 'sphere')
            single_node_idx = 744;
        else
            single_node_idx = 29584;
        end
        sn.Qe(k,:) = last.Qe(single_node_idx);
        sn.Ve(k,:) = last.Ve(single_node_idx);
        sn.dVe(k,:) = last.dVe(single_node_idx);
        sn.Qi(k,:) = last.Qi(single_node_idx);
        sn.Vi(k,:) = last.Vi(single_node_idx);
        sn.dVi(k,:) = last.dVi(single_node_idx);
        sn.Qi_fs(k,:) = last.Qi_fs(single_node_idx);
        sn.Vi_fs(k,:) = last.Vi_fs(single_node_idx);
        sn.dVi_fs(k,:) = last.dVi_fs(single_node_idx);
        sn.D22(k,:) = last.D22(single_node_idx);
        sn.K(k,:) = last.K(single_node_idx);
    end

    % plot frame for video and write frame
    if flag_video
        clf(f);
        if strcmp(type, 'sphere')
            plot_sphere_instance(locs, last, macro_pos, micro_pos);
        else
            plot_brain_instance(surf, surf_sphere, last, macro_pos, micro_pos);
        end
        drawnow;
        im = getframe(f);
        writeVideo(vidObj,im);
    end
end

save(SAMPLE_DATA_FILE, ...
    'Qe_rand', 'Ve_rand', 'Qe_avg', 'Ve_avg', ...
    'Qe_macro', 'Ve_macro', 'Qe_micro', 'Ve_micro', 'sn');

if flag_video
    close(vidObj);
end
% allocate storing space
Qe_rand = NaN(K*T, 3);
Ve_rand = NaN(K*T, 3);
Qe_avg = NaN(K*T, 3);
Ve_avg = NaN(K*T, 3);
Qe_macro = NaN(K*T, 19);
Ve_macro = NaN(K*T, 19);
Qe_micro = NaN(K*T, 19);
Ve_micro = NaN(K*T, 19);

% load node positions
load('N10242_R10_wideNodes.mat');
pos_hemi = locs(:,3) >= 0;
neg_hemi = locs(:,3) < 0;

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

    Qe_rand(1+(k-1)*T : k*T,1) = fine.Qe_focus(:,1);
    Qe_rand(1+(k-1)*T : k*T,2) = fine.Qe_macro(:,1);
    Qe_rand(1+(k-1)*T : k*T,3) = fine.Qe_macro(:,8);
    Qe_rand(1+(k-1)*T : k*T,4) = fine.Qe_normal(:,1);
    Ve_rand(1+(k-1)*T : k*T,1) = fine.Ve_focus(:,1);
    Ve_rand(1+(k-1)*T : k*T,2) = fine.Ve_macro(:,1);
    Ve_rand(1+(k-1)*T : k*T,3) = fine.Ve_macro(:,8);
    Ve_rand(1+(k-1)*T : k*T,4) = fine.Ve_normal(:,1);

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
    

    % plot frame for video and write frame
    scatter(locs(pos_hemi,1), locs(pos_hemi,2), 15, last.Qe(pos_hemi), 'filled');
    caxis([0,30]);
    xlim([-10,33]);
    
    hold on;
    scatter(23 - locs(neg_hemi,1), locs(neg_hemi,2), 15, last.Qe(neg_hemi), 'filled');
    caxis([0,30]);
    drawnow;
        
    scatter(23 - locs(micro_idx,1),locs(micro_idx,2), 15, 'g', 'filled');
    scatter(23 - locs(macro_idx,1), locs(macro_idx,2), 15, 'r', 'filled');
    
    f = getframe;
    writeVideo(vidObj,f);
    clf;
end

save(SAMPLE_DATA_FILE, 'Qe_rand', 'Ve_rand', 'Qe_avg', 'Ve_avg', ...
    'Qe_macro', 'Ve_macro', 'Qe_micro', 'Ve_micro');

close(vidObj);
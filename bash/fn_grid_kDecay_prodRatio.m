function fn_grid_kDecay_prodRatio(kDecay, prodRatio)

    DATA_ROOT_DIR = './data/grid_kDecay_prodRatio/';
    mkdir(DATA_ROOT_DIR);

    %% load grid
    load('N10242_R10.mat');
    
    %% initialize parameters and map
    K = 3000;
    T0 = 0.1;
    map = make_map(laplacian);

    %% initialize initial state
    last = make_IC(N);
    last.Qi_fs = last.Qi;
    last.Vi_fs = last.Vi;
    last.F_ii_fs = last.F_ii;
    last.Phi_ii_fs = last.Phi_ii;
    last.dVi_fs = last.dVi;

    %% define zones
    lessihb_filter = lessihb_area;

    zones.focus_zone = map == 1;
    zones.lessihb_zone = lessihb_filter & map ~= 1;
    zones.normal_zone = ~lessihb_filter;

    normal_sample_idx = [];
    macro_idx = [744 437 821 1141 1140 820 436 251 555 981 1253 1585 1537 1584 1252 980 554 250 187];
    micro_idx = [744 659 753 837 836 752 658 579 669 777 845 933 929 932 844 776 668 578 573];
    fine_idx = union(find(map), [macro_idx, micro_idx]);

    %% initialize constants and make modifications
    global HL
    HL = SCM_init_globs(N);

    HL.kR = 7 * ones(N,1);
    HL.kR(zones.normal_zone) = 0;

    HL.prodRatio = prodRatio;
    HL.k_decay = kDecay;

    HL.KtoVe = 0;
    HL.KtoVi = 0;
    HL.KtoVi_fs = 0;
    HL.KtoD  = -2;
    HL.D22min = 0.1;
    HL.FS_ratio = 0;

    last.D22(:) = 7; last.D11 = last.D22/100;
    last.K(:) = 5;
    last.K(map == 1) = 12;
    
    %% set the output directory and save meta file
    id = ['kDecay' num2str(kDecay, '%.1f') '_prodRatio' num2str(prodRatio, '%.1f')];
    folder_name = ['sphere_N' num2str(N) '_R' num2str(R) '_' id];

    OUTPUT_DIR = [DATA_ROOT_DIR folder_name '/raw/'];
    mkdir(OUTPUT_DIR);
    
    META_FILE = [DATA_ROOT_DIR folder_name '/vars.mat'];
    save(META_FILE, 'HL', 'map', 'fine_idx', 'normal_sample_idx', 'last');

    %% run simulation
    for k = 1:K
        source_drive = NaN;
        tic;

        [samp_time,last,fine] = seizing_cortical_field(...
            source_drive, map, T0, last, ...
            locs, laplacian, avg_D, ...
            zones, fine_idx, normal_sample_idx, ...
            true);

        save([OUTPUT_DIR 'seizing_cortical_field_k_' num2str(k) '.mat'], ...
            'samp_time', 'last', 'fine');
    end

    %% run analysis
    [~,macro_speed,micro_speed,recruitment_speed] = ...
        main_plot_graphs(id, DATA_ROOT_DIR, true, false, true);
    
    SPEED_FILE = [DATA_ROOT_DIR folder_name '/speeds.mat'];
    save(SPEED_FILE, 'macro_speed', 'micro_speed', 'recruitment_speed');
    
    %% print speeds to a file
    fileID = fopen([DATA_ROOT_DIR 'output.txt'],'a');
    fprintf(fileID, '%s --- %.3f %.3f %.3f\n', ...
        id, macro_speed, micro_speed, recruitment_speed);
    fclose(fileID);
end
function fn_grid_FStwo_source(source_arg)

    DATA_ROOT_DIR = './data/grid_FStwo_source/';
    mkdir(DATA_ROOT_DIR);

    %% load grid
    load('N10242_R10.mat');
    
    %% initialize parameters and map
    K = 2000;
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
    lessihb_filter = true(N, 1);

    zones.focus_zone = map == 1;
    zones.lessihb_zone = lessihb_filter & map ~= 1;
    zones.normal_zone = ~lessihb_filter;

    normal_sample_idx = 10242;
    
    macro_idx = [744 437 821 1141 1140 820 436 251 555 981 1253 1585 1537 1584 1252 980 554 250 187];
    micro_idx = [744 659 753 837 836 752 658 579 669 777 845 933 929 932 844 776 668 578 573];
    fine_idx = union(find(map), [macro_idx, micro_idx]);

    %% initialize constants and make modifications
    global HL
    HL = SCM_init_globs(N);

    HL.kR = 2.5;
    HL.KtoVe = 0;
    HL.KtoVi = 0;
    HL.KtoVi_fs = 2000;
    HL.KtoD  = -2;
    HL.D22min = 0.1;

    last.D22(:) = 6; last.D11 = last.D22/100;

    %% set the output directory and save meta file
    id = ['source' num2str(source_arg)];
    folder_name = ['sphere_N' num2str(N) '_R' num2str(R) '_' id];

    OUTPUT_DIR = [DATA_ROOT_DIR folder_name '/raw/'];
    mkdir(OUTPUT_DIR);

    META_FILE = [DATA_ROOT_DIR folder_name '/vars.mat'];
    save(META_FILE, 'HL', 'map', 'lessihb_idx', 'normal_sample_idx', 'last');

    %% run simulation
    for k = 1:K
        source_drive = source_arg;
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
    main_plot_graphs(id, DATA_ROOT_DIR, true, false, true);
end
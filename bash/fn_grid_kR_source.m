function fn_grid_kR_source(kR_arg, source_arg, D22min_arg)

    DATA_ROOT_DIR = './data/grid_kR_source/';
    mkdir(DATA_ROOT_DIR);

    %% load grid
    load('N10242_R10.mat');
    
    %% initialize parameters and map
    K = 200;
    T0 = 1;
    map = make_map(laplacian);

    %% initialize initial state
    last = make_IC(N);

    %% define zones
    lessihb_filter = locs(:,3) < -6;

    zones.focus_zone = map == 1;
    zones.lessihb_zone = lessihb_filter & map ~= 1;
    zones.normal_zone = ~lessihb_filter;

    lessihb_idx = find(lessihb_filter);
    normal_sample_idx = 500;

    %% initialize constants and make modifications
    global HL
    HL = SCM_init_globs(N);

    HL.kR = kR_arg;
    HL.KtoVe = 0;
    HL.KtoVi = 0;
    HL.KtoD  = -20;
    HL.D22min = D22min_arg;

    %% set the output directory and save meta file
    id = ['kR' num2str(kR_arg) '_source' num2str(source_arg) '_D22min' num2str(D22min_arg)];
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
            zones, lessihb_idx, normal_sample_idx, ...
            true);

        save([OUTPUT_DIR 'seizing_cortical_field_k_' num2str(k) '.mat'], ...
            'samp_time', 'last', 'fine');
    end

    %% run analysis
    main_plot_graphs(id, DATA_ROOT_DIR, false, false);
end
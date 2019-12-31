function [fg_joint,macro_speed,micro_speed,recruitment_speed] = ...
    main_plot_graphs(id, DATA_ROOT_DIR, flag_plot, flag_video, flag_tabbed)

    %% figure properties
    set(0, 'DefaultTextFontsize', 8);
    set(0, 'DefaultFigurePosition', [600, 50, 1000, 900]);

    % brain type only
    NOTE = 'closest7_dipoles_window15s_r1dist';
    % loc_grid_center = [67.83, -28.99, 27.62];
    loc_grid_center = [64.41, -7.28, 21.48];
    dist_grid_macro = 9;                                  % distance between macro-electrodes (mm)
    dist_grid_micro = 0.5;                                 % distance between micro-electrodes (mm)
    flag_dipole = true;
    closest_N = 7;

    %% *** SPECIFY *** input directory and file names
    [DATA_DIR, type] = find_full_id(DATA_ROOT_DIR, id);
    RAW_DIR = [DATA_DIR 'raw/'];
    META_FILE = [DATA_DIR 'vars.mat'];
    load(META_FILE);

    %% *** SPECIFY *** output directory and file names
    if strcmp(type, 'sphere')
        ANALYSIS_DIR = DATA_DIR;
    else
        ANALYSIS_DIR = [DATA_DIR 'grid_'...
                                 num2str(loc_grid_center(1), '%.2f') '_' ...
                                 num2str(loc_grid_center(2), '%.2f') '_' ...
                                 num2str(loc_grid_center(3), '%.2f') '_' ...
                                 num2str(dist_grid_macro) '_' num2str(dist_grid_micro) '_' NOTE '/'];
        mkdir(ANALYSIS_DIR);
    end

    VIDEO_FILE = [ANALYSIS_DIR 'movie.mp4'];
    SAMPLE_DATA_FILE = [ANALYSIS_DIR 'data_sample.mat'];
    COHERENCE_FILE   = [ANALYSIS_DIR 'data_coherence_dense.mat'];
    ELEC_FILE        = [ANALYSIS_DIR 'data_electrode.mat'];

    %% *** SPECIFY *** figure names
    if flag_tabbed && flag_plot
        JOINT_FIG    = [ANALYSIS_DIR 'joint.fig'];
        fg_joint = figure('Name', 'Joint');
        tgroup = uitabgroup(fg_joint);

        tab_course_neg = uitab(tgroup, 'Title', 'Course (neg)');
        tab_course_pos = uitab(tgroup, 'Title', 'Course (pos)');
        tab_traces = uitab(tgroup, 'Title', 'Traces');
        tab_single = uitab(tgroup, 'Title', 'Single Node');
        tab_coh_macro = uitab(tgroup, 'Title', 'Coherence (macro)');
        tab_coh_micro = uitab(tgroup, 'Title', 'Coherence (micro)');

        [fg_course_neg, fg_course_pos, fg_traces, ...
            fg_single, fg_coh_macro, fg_coh_micro] = deal(fg_joint);
    
    elseif flag_plot
        COURSE_NEG_FIG   = [ANALYSIS_DIR 'fig_course_neg.fig'];
        COURSE_POS_FIG   = [ANALYSIS_DIR 'fig_course_pos.fig'];
        TRACES_FIG       = [ANALYSIS_DIR 'fig_traces.fig'];
        SINGLE_FIG       = [ANALYSIS_DIR 'fig_single.fig'];
        COH_MACRO_FIG    = [ANALYSIS_DIR 'fig_coh_summary_macro.fig'];
        COH_MICRO_FIG    = [ANALYSIS_DIR 'fig_coh_summary_micro.fig'];

        fg_joint = NaN;
        fg_course_neg = figure('Name', 'Course (neg)'); tab_course_neg = gcf;
        fg_course_pos = figure('Name', 'Course (pos)'); tab_course_pos = gcf;
        fg_traces = figure('Name', 'Traces'); tab_traces = gcf;
        fg_single = figure('Name', 'Single Node'); tab_single = gcf;
        fg_coh_macro = figure('Name', 'Coherence (macro)'); tab_coh_macro = gcf;
        fg_coh_micro = figure('Name', 'Coherence (micro)'); tab_coh_micro = gcf;
    end

    %% *** SPECIFY *** time sequences
    K = length(dir([RAW_DIR 'seizing_*.mat']));
    T0 = 1 * (K < 1000) + 0.1 * (K >= 1000);
    T = 500 * T0;

    total_time = K*T0;
    sparse_time = (1:K) * T0;
    fine_time = (1:K*T) * (T0/T);

    % for coherence split the entire course into P periods
    wind = 10; % seconds
    P = total_time - (wind - 1);
    per_P = wind * 500;

    %% *** SPECIFY *** grid
    if strcmp(type, 'sphere')
        load('N10242_R10.mat');
        pos_hemi = locs(:,3) >= 0; %#ok<*NODEF,*NASGU>
        neg_hemi = locs(:,3) < 0;
    else
        load('N40962.mat');
        load('unitsphere.mat');

        surf.vertices = locs;
        surf.faces = tri;
        surf_sphere.vertices = 10 * coord';
        surf_sphere.faces = tri; %#ok<*STRNU>
    end

    %% *** LOAD/GENERATE *** sample data
    if exist(SAMPLE_DATA_FILE, 'file') == 2 && exist(ELEC_FILE, 'file') == 2
        load(ELEC_FILE);
        load(SAMPLE_DATA_FILE);
    else
        prepare_electrode_data;
        prepare_sample_data;
    end

    %% *** LOAD/GENERATE *** coherence data
    if exist(COHERENCE_FILE, 'file') == 2
        load(COHERENCE_FILE);
    else
        prepare_coherence;
    end

    if flag_plot
        %% *** PLOT *** seizure course
        plot_course;

        %% *** PLOT *** firing rate and voltage traces
        plot_traces;

        %% *** PLOT *** single node dynamics
        figure(fg_single); axes('parent', tab_single);
        
        % legend('Qe', 'Qi', 'Ve', 'Vi', 'D22', 'dVe', 'dVi', 'K', ...
        %     'Qi\_fs', 'Vi\_fs', 'dVi\_fs');
        
        subplot(2,2,1);
        plot(sparse_time, single_node(:,[2,1,9,8]));
        legend('Qi', 'Qe', 'Qi\_fs', 'K');
        xlim([25,300]);
        xlabel('time (s)');
        ylabel('firing rate (Hz)');
        
        subplot(2,2,2);
        plot(sparse_time, single_node(:,[2,1,9,8]));
        legend('Qi', 'Qe', 'Qi\_fs', 'K');
        xlim([270,300]);
        xlabel('time (s)');
        ylabel('firing rate (Hz)');

        %% *** PLOT *** coherence statistics
        central_t = (wind/2) : (wind/2 + P - 1);
        [~,period_idx] = min(abs(central_t - 240)); % estimate wave for this period

        warning('off', 'stats:statrobustfit:IterationLimit');
        
        [t_coh,t_coh_conf,t_phi,electrode_2d] = deal( ...
            macro_t_coh, macro_t_coh_conf, macro_t_phi, macro_2d);
        fprintf('macro electrodes ');
        figure(fg_coh_macro); axes('parent', tab_coh_macro);
        plot_coherence;
        macro_speed = swd_speed;

        if ~isempty(micro_pos)
            [t_coh,t_coh_conf,t_phi,electrode_2d] = deal( ...
                micro_t_coh, micro_t_coh_conf, micro_t_phi, micro_2d); %#ok<*ASGLU>
            fprintf('micro electrodes ');
            figure(fg_coh_micro); axes('parent', tab_coh_micro);
            plot_coherence;
            micro_speed = swd_speed;
        else
            micro_speed = NaN;
        end

        %% *** CALCULATE *** recruitment speed
        
        node_to_focus = dist([locs(focus_idx(1),:); macro_pos]');
        node_to_focus = node_to_focus(1,2:end);
        
        [~, node_e] = min(node_to_focus);
        [~, node_l] = max(node_to_focus);

        t_e = fine_time(find(Qe_macro(:,node_e) > 15, 1));
        t_l = fine_time(find(Qe_macro(:,node_l) > 15, 1));

        travel_dist = dist(macro_2d([node_e, node_l], :)');
        if ~isempty(t_e) && ~isempty(t_l)
            recruitment_speed = travel_dist(1,2) / (t_l-t_e);
        else
            recruitment_speed = NaN;
        end
        fprintf(['recruitment speed is ' num2str(recruitment_speed) ' cm/s \n']);

        %% *** SAVE *** figures
        if flag_tabbed
            saveas(fg_joint, JOINT_FIG);
        else
            saveas(fg_course_neg, COURSE_NEG_FIG);
            saveas(fg_course_pos, COURSE_POS_FIG);
            saveas(fg_traces, TRACES_FIG);
            saveas(fg_single, SINGLE_FIG);
            saveas(fg_coh_macro, COH_MACRO_FIG);
            saveas(fg_coh_micro, COH_MICRO_FIG);
        end
    end
end
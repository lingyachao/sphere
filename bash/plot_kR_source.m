DATA_ROOT_DIR = './data/grid_kR_source/';

kR_arg_list = [1.5 2.5 3.5];
source_arg_list = [2 2.5 3 3.5];
D22min_arg_list = [0.1 0.6 1.1];

[grid_macro_speed, grid_micro_speed, grid_recruitment_speed] = ...
    deal(NaN(length(kR_arg_list), length(source_arg_list), length(D22min_arg_list)));

for kR_idx = 1:length(kR_arg_list)
    for source_idx = 1:length(source_arg_list)
        for D22min_idx = 1:length(D22min_arg_list)
            
            kR_arg = kR_arg_list(kR_idx);
            source_arg = source_arg_list(source_idx);
            D22min_arg = D22min_arg_list(D22min_idx);
            
            id = ['kR' num2str(kR_arg) '_source' num2str(source_arg) '_D22min' num2str(D22min_arg)];
            fprintf(['analyzing...' id '\n']);
            
            [fg_joint, ...
                macro_speed_grid(kR_idx, source_idx, D22min_idx), ...
                micro_speed_grid(kR_idx, source_idx, D22min_idx), ...
                recruitment_speed(kR_idx, source_idx, D22min_idx)] = ...
                main_plot_graphs(id, DATA_ROOT_DIR, true, true); %#ok<*SAGROW>
            
            saveas(fg_joint, [DATA_ROOT_DIR 'joint_figures/' id '.fig']);
            save([DATA_ROOT_DIR 'grid_speeds.mat'], ...
                'grid_macro_speed', 'grid_micro_speed', 'grid_recruitment_speed');
        end
    end
end
if strcmp(type, 'sphere')
    [focus_idx, macro_pos, macro_transform, macro_2d, ...
                micro_pos, micro_transform, micro_2d] = ...
        generate_electrode_grid_sphere(RAW_DIR);
else
    % find_electrode_grid_center(RAW_DIR, dist_grid);
    % keyboard;
    
    [focus_idx, macro_pos, macro_transform, macro_2d, ...
                micro_pos, micro_transform, micro_2d] = ...
        generate_electrode_grid_brain(loc_grid_center, dist_grid_macro, dist_grid_micro, RAW_DIR, ...
                                      flag_dipole, closest_N);
end

save(ELEC_FILE, 'focus_idx', 'macro_pos', 'macro_transform', 'macro_2d', ...
                             'micro_pos', 'micro_transform', 'micro_2d');
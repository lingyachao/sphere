function [focus_idx, macro_pos, macro_transform, macro_2d, ...
                     micro_pos, micro_transform, micro_2d] = generate_electrode_grid_sphere(RAW_DIR)

    load('N10242_R10.mat');
    load('node_2d_N19.mat');
    
    focus_idx = find(make_map(laplacian));
    
    macro_idx = [744 437 821 1141 1140 820 436 251 555 981 1253 1585 1537 1584 1252 980 554 250 187];
    macro_transform = zeros(19, N);
    macro_transform(sub2ind(size(macro_transform), 1:19, macro_idx)) = 1;
    macro_pos = locs(macro_idx,:);
    macro_2d = position;
    
    micro_idx = [744 659 753 837 836 752 658 579 669 777 845 933 929 932 844 776 668 578 573];
    micro_transform = zeros(19, N);
    micro_transform(sub2ind(size(micro_transform), 1:19, micro_idx)) = 1;
    micro_pos = locs(micro_idx,:);
    micro_2d = position;
    
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(50) '.mat']);
    f = figure;
    set(f, 'Position', [200 300 900 400]);
    plot_sphere_instance(locs, last, macro_pos, micro_pos);
end











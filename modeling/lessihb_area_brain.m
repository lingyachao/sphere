function lessihb_idx = lessihb_area_brain

    load('./computed_sphere_grid/unitsphere.mat', 'coord');
    
    boundary_idx = [123 985 1028 2025 2401 2402 2433 2435 4010 4043 4049 4167 4173 4207 6959 6985 7079 7099 7182 7199 7200 7391 7680 7684 7730 7812 7847 7866 8051 9574 9579 9582 9692 9735 9749 9913 9944 9955 15357 15359 15394 15395 15407 15408 15920 16048 16083 16092 16097 16098 16099 16241 16242 16246 16559 16561 16564 16566 16589 16591 16597 16714 16721 16725 27720 27724 27824 27831 27832 27851 27852 27897 27900 27909 27916 28117 28194 28206 28209 28214 28227 28229 28272 28273 28426 28606 28608 28618 28620 28643 28646 29266 29267 29269 29394 29401 29404 29423 29427 29435 29452 30528 30541 30544 30575 30577 30596 30606 30618 30621 30650 30793 30801 30805 30914 30916 30933 30934 30936 30941 30942 31263 31269 31271 31314 31316 31323 31327 31344 31349 31360 31938 31942 31945 31959 31998 32020 32023 32051 32115 32128 38059 38156 38157 38159 38208 38231 38241 38434 38449 38452 38556 38670 38672 38679 38681 38695 38709 38846 38847 38849 38898 38903 38907 38919 39065 39593 39597 39619 39622 39626 39627 39628 39642 39721 39733 39736 39764 39771];
    k = boundary(coord(2,boundary_idx)', coord(3,boundary_idx)', 0.8);
    boundary_idx = boundary_idx(k);
    
    p = coord(1,:)' > 0;
    
    % scatter(coord(2,:)', coord(3,:)', 15, 'y', 'filled');
    % hold on;
    % plot(coord(2,boundary_idx)', coord(3,boundary_idx)');
    
    [in, on] = inpolygon(coord(2,:)', coord(3,:)', ...
        coord(2,boundary_idx)', coord(3,boundary_idx)');
    lessihb_idx = (in | on) & p;
    % scatter(coord(2,lessihb_idx)', coord(3,lessihb_idx)', 15, 'g', 'filled');
    
end

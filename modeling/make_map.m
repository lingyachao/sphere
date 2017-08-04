function map = make_map(laplacian)
% make a map for the source

    map = 1 * (laplacian(:,15473) ~= 0);
end
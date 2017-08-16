function map = make_map(laplacian)
% make a map for the source
    if size(laplacian, 1) > 20000
        map = 1 * (laplacian(:,15473) ~= 0);
    else
        map = 1 * (laplacian(:,1) ~= 0);
    end
end
function map = make_map(laplacian)
% make a map for the source
    if size(laplacian, 1) > 20000
        map = 1 * (laplacian(:,35) ~= 0);
        % map = 1.2 * (laplacian(:,15473) ~= 0) + 1 * (laplacian(:,30321) ~= 0);
    else
        map = 1 * (laplacian(:,1) ~= 0);
    end
end
function wl = get_weak_laplacian(locs, laplacian)

    wl = laplacian;

    z_thres = -8;
    strength = 0.1;
    
    wl(locs(:,3) > z_thres, locs(:,3) < z_thres) = ...
        strength * wl(locs(:,3) > z_thres, locs(:,3) < z_thres);
    wl(locs(:,3) < z_thres, locs(:,3) > z_thres) = ...
        strength * wl(locs(:,3) < z_thres, locs(:,3) > z_thres);

    wl = sparse(wl + (-1)*diag(sum(wl)));
end
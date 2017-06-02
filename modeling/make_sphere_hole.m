function [N,locs,laplacian,avg_D] = make_sphere_hole(targetN, R)
% create a sphere grid with number of nodes closest to targetN
% output includes the cartesian locations, the average distance between
% nodes and the Laplacian matrix operator

    [lat,long] = GridSphere(targetN);
    [x,y,z] = sph2cart(deg2rad(long), deg2rad(lat), R);
    locs = [x,y,z];
    N = length(lat);

    D = squareform(pdist([x,y,z]));
    D = D + max(max(D)) * eye(size(D));

    laplacian = zeros(N, N);
    neigh_D = NaN(N, 6);
    
    for i = 1 : N
        [sorted,I] = sort(D(i,:));
        I = I(1:6);
        laplacian(i,I) = 1;
        neigh_D(i,:) = sorted(1:6);
    end

    flt = z < 0.9 * R;
    N = sum(flt);
    locs = locs(flt,:);
    laplacian = laplacian(flt,flt);
    
    avg_D = mean(mean(neigh_D));
    laplacian = sparse(laplacian + (-1)*diag(sum(laplacian)));
    % laplacian = sparse(laplacian + (-6)*eye(N));
    
%     deg = 5;
%     micro_lats = [deg, deg, deg, 0, 0, 0, -deg, -deg, -deg]';
%     micro_longs = [-deg, 0, deg, -deg, 0, deg, -deg, 0, deg]';
%     
%     micro_idx = find_neighbors(micro_lats, micro_longs, lat, long);
%     macro_idx = find_neighbors(micro_lats * 4, micro_longs * 4, lat, long);
end

% function idx = find_neighbors(q_lats, q_longs, all_lats, all_longs)
% 
%     idx = NaN(length(q_lats), 1);
%     
%     for i = 1 : length(q_lats)
%         arcs = NaN(length(all_lats), 1);
%         for j = 1 : length(all_lats)
%             arcs(j) = distance(q_lats(i), q_longs(i), all_lats(j), all_longs(j));
%         end
%         [~,I] = min(arcs);
%         idx(i) = I;
%     end
% 
% end

    


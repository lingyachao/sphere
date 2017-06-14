function [N,locs,laplacian,avg_D] = make_brain_grid

    load('autism.surface.mat', 'controlouter', 'tri');
    
    N = size(controlouter, 3);
    locs = NaN(N, 3);
    locs(:,1) = controlouter(1,1,:);
    locs(:,2) = controlouter(1,2,:);
    locs(:,3) = controlouter(1,3,:);
    
    laplacian = sparse(zeros(N, N));
    
    for i = 1 : length(tri)
        laplacian(tri(i,1), tri(i,[2,3])) = 1;
        laplacian(tri(i,2), tri(i,[1,3])) = 1;
        laplacian(tri(i,3), tri(i,[1,2])) = 1;
    end

    avg_D = 0.3;
    laplacian = sparse(laplacian + (-1)*diag(sum(laplacian)));
end
% create operator for zero-flux
% edge = sum(laplacian, 2) ~= 0;
% avg_flux = laplacian(edge,:);
% avg_flux = max(avg_flux, 0);
% avg_flux = avg_flux ./ repmat(sum(avg_flux, 2), 1, N);
% avg_flux = sparse(zeros(sum(edge), N));
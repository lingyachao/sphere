[DATA_DIR, type] = find_full_id('./data/', '09211929');
RAW_DIR = [DATA_DIR 'raw/'];
K = length(dir([RAW_DIR 'seizing_*.mat']));
N = 10242;

load('N10242_R10.mat');

Qe_all = NaN(N, K);
recruited = NaN(N, 1);

for k = 1:K
    fprintf(['Read in ' num2str(k) '\n']);
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(k) '.mat'], 'last');
    Qe_all(:,k) = last.Qe;
end

for n = 1:N
    recruited_t = find(Qe_all(n,:) > 15, 1);
    if ~isempty(recruited_t)
        recruited(n) = recruited_t;
    end
end

[lat,long] = GridSphere(10242);
arc_dist = 10 * (lat + 90) * pi/180;
plot(arc_dist, recruited);

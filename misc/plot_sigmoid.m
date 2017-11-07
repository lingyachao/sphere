
figure;
Vi_grid = -80:0.1:0;
Qi_grid = HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i/2) .* (Vi_grid - HL.theta_i)))) ...     % The I voltage must be big enough,
                  - HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i/2) .* (Vi_grid - (HL.theta_i + 15)))));     % ... but not too big.
plot(Vi_grid, Qi_grid, 'k', 'LineWidth', 1);
xlabel('Vi');
ylabel('Qe/Qi/Qi_{fs}');

hold on;

Qi_grid = HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i) .* (Vi_grid - HL.theta_i)))) ...     % The I voltage must be big enough,
                  - HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i) .* (Vi_grid - (HL.theta_i+25)))));     % ... but not too big.
plot(Vi_grid, Qi_grid, 'k', 'LineWidth', 1);



figure;
K = 0:0.1:15;
PK = 0.5;
ghk = 26.7123 * log((PK*K + 11.75) ./ (PK*(150-10*K) + 50.25));
plot(K, ghk)
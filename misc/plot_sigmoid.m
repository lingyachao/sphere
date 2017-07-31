
Vi_grid = -80:0;
Qi_grid = HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i) .* (Vi_grid - HL.theta_i)))) ...     % The I voltage must be big enough,
                  - HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i) .* (Vi_grid - (HL.theta_i+30)))));     % ... but not too big.
plot(Vi_grid, Qi_grid);              
function [samp_time,last,fine] = seizing_cortical_field( ...
    source_del_VeRest, map, time_end, IC, ...
    locs, lap_synap, lap_diff, avg_D, ...
    zones, fine_idx, normal_sample_idx, ...
    save_output)
 
    global HL
    
    % set no. of sampling points (must be even!) along each axis of cortical grid
    N = size(lap_diff, 1);
    
    % initialize random number generator (from input argument)
    rand_state = sum(100*clock);
    rng(rand_state); % randn('state', rand_state);

    % set time resolution
    dt = 0.2 * 1e-3;

    % number of time-steps for simulation
    Nsteps = round(time_end/dt);

    % use as initial conditions the "last" values of previous simulation.
    Qe_grid = IC.Qe;
    Qi_grid = IC.Qi;
    Qi_fs_grid = IC.Qi_fs;
    Ve_grid = IC.Ve;
    Vi_grid = IC.Vi;
    Vi_fs_grid = IC.Vi_fs;

    phi_ee = IC.phi_ee;
    phi_ei = IC.phi_ei;
    phi2_ee = IC.phi2_ee;
    phi2_ei = IC.phi2_ei;

    F_ee = IC.F_ee;
    F_ei = IC.F_ei;
    F_ie = IC.F_ie;
    F_ii = IC.F_ii;
    F_ii_fs = IC.F_ii_fs;

    Phi_ee = IC.Phi_ee;
    Phi_ei = IC.Phi_ei;
    Phi_ie = IC.Phi_ie;
    Phi_ii = IC.Phi_ii;
    Phi_ii_fs = IC.Phi_ii_fs;

    D11 = IC.D11;
    D22 = IC.D22;
    del_VeRest = IC.dVe;
    del_ViRest = IC.dVi;
    del_ViRest_fs = IC.dVi_fs;
    K = IC.K;
        
    samp_rate = 10;
    samp_time = (0:samp_rate:Nsteps-1)' * dt;
    N_samp = Nsteps/samp_rate;

    if save_output
        % allocate matrices for recording Qe/Ve (fine time scale) at micro/macro
        % nodes, focus nodes, normal nodes
    
        fine.Qe_lessihb = zeros(N_samp, length(fine_idx));
        fine.Qe_normal = zeros(N_samp, length(normal_sample_idx));

        fine.Qe_focus_avg = zeros(N_samp, 1);
        fine.Qe_lessihb_avg = zeros(N_samp, 1);
        fine.Qe_normal_avg = zeros(N_samp, 1);

        fine.Ve_lessihb = zeros(N_samp, length(fine_idx));
        fine.Ve_normal = zeros(N_samp, length(normal_sample_idx));

        fine.Ve_focus_avg = zeros(N_samp, 1);
        fine.Ve_lessihb_avg = zeros(N_samp, 1);
        fine.Ve_normal_avg = zeros(N_samp, 1);
    
    else
        fine = NaN;
    end
    
    % noise-amplitude coefficients for subcortical flux (note 1/sqrt(dt) factor)
    B_ee = HL.noise_sf * sqrt(HL.noise_sc * HL.phi_ee_sc / dt);
    B_ei = HL.noise_sf * sqrt(HL.noise_sc * HL.phi_ei_sc / dt);
    
    for i = 0:Nsteps-1

        % 0. record Qe/Ve every 10 steps
        if save_output && mod(i, samp_rate) == 0
            idx = i/samp_rate + 1;
            
            fine.Qe_lessihb(idx,:) = Qe_grid(fine_idx);
            fine.Qe_normal(idx,:) = Qe_grid(normal_sample_idx);

            fine.Qe_focus_avg(idx,:) = mean(Qe_grid(zones.focus_zone));
            fine.Qe_lessihb_avg(idx,:) = mean(Qe_grid(zones.lessihb_zone));
            fine.Qe_normal_avg(idx,:) = mean(Qe_grid(zones.normal_zone));

            fine.Ve_lessihb(idx,:) = Ve_grid(fine_idx);
            fine.Ve_normal(idx,:) = Ve_grid(normal_sample_idx);

            fine.Ve_focus_avg(idx,:) = mean(Ve_grid(zones.focus_zone));
            fine.Ve_lessihb_avg(idx,:) = mean(Ve_grid(zones.lessihb_zone));
            fine.Ve_normal_avg(idx,:) = mean(Ve_grid(zones.normal_zone));
        end
        
        % 1. update wave equations
        phi2_ee_1 = phi2_ee + dt * (-2 * HL.v * HL.Lambda * phi2_ee ...
                            - (HL.v * HL.Lambda)^2 * phi_ee ...
                            + (HL.v * HL.Lambda)^2 * Qe_grid)...
                            + dt * (HL.v ./ avg_D) .^2 .* (lap_synap * phi_ee);
        phi_ee_1  = phi_ee + dt * phi2_ee;

        phi2_ei_1 = phi2_ei + dt * (-2 * HL.v * HL.Lambda * phi2_ei ...
                            - (HL.v * HL.Lambda)^2 * phi_ei ...
                            + (HL.v * HL.Lambda)^2 * Qe_grid) ...
                            + dt * (HL.v ./ avg_D) .^2 .* (lap_synap * phi_ei);
        phi_ei_1  = phi_ei + dt * phi2_ei;

        % 2. update the 4 synaptic flux equations (include sc noise)

        %%%% E-to-E %%%%
        F_ee_1   = F_ee + dt * HL.gamma_e^2 * (-2/HL.gamma_e*F_ee - Phi_ee ...
                        + HL.Nee_a * phi_ee ...             % long range
                        + HL.Nee_b * Qe_grid ...            % short range
                        + HL.noise_sc * HL.phi_ee_sc ...    % subcortical (tonic)
                        + B_ee .* randn(N, 1));              % subcortical (random)
        Phi_ee_1 = Phi_ee + dt*F_ee;

        %%%% E-to-I %%%%
        F_ei_1   = F_ei + dt * HL.gamma_e^2 * (-2/HL.gamma_e*F_ei - Phi_ei ...
                        + HL.Nei_a * phi_ei ...             % long range
                        + HL.Nei_b * Qe_grid ...            % short range
                        + HL.noise_sc * HL.phi_ei_sc ...    % subcortical (tonic)
                        + B_ei * randn(N, 1));              % subcortical (random)
        Phi_ei_1 = Phi_ei + dt*F_ei;

        %%%% I-to-E %%%%
        F_ie_1   = F_ie + dt * HL.gamma_i^2 * (-2/HL.gamma_i*F_ie - Phi_ie ...
                        + HL.Nie_b .* Qi_grid ...
                        + HL.Nie_fs .* Qi_fs_grid); % short range
        Phi_ie_1 = Phi_ie + dt*F_ie;

        %%%% I-to-I %%%%
        F_ii_1   = F_ii + dt * HL.gamma_i^2 * (-2/HL.gamma_i*F_ii - Phi_ii ...
                        + HL.Nii_b .* Qi_grid); % short range
        Phi_ii_1 = Phi_ii + dt*F_ii;
        
        F_ii_fs_1   = F_ii_fs + dt * HL.gamma_i^2 * (-2/HL.gamma_i*F_ii_fs - Phi_ii_fs ...
                        + HL.Nii_fs .* Qi_fs_grid); % short range
        Phi_ii_fs_1 = Phi_ii_fs + dt*F_ii_fs;

        % 3. update the soma voltages

        Ve_grid_1 = Ve_grid + dt/HL.tau_e * ((HL.Ve_rest - Ve_grid) + del_VeRest ...
                            + HL.ge .* Psi_ee(Ve_grid) .* Phi_ee ...      %E-to-E
                            + HL.gi * Psi_ie(Ve_grid) .* Phi_ie ...      %I-to-E
                            + D11 .* (lap_synap * Ve_grid));

        Vi_grid_1 = Vi_grid + dt/HL.tau_i * ((HL.Vi_rest - Vi_grid) + del_ViRest ...
                            + HL.ge .* Psi_ei(Vi_grid) .* Phi_ei ...      %E-to-I
                            + HL.gi * Psi_ii(Vi_grid) .* Phi_ii ...      %I-to-I
                            + D22 .* (lap_synap * Vi_grid));
        
        Vi_fs_grid_1 = Vi_fs_grid + dt/HL.tau_i * ((HL.Vi_rest_fs - Vi_fs_grid) + del_ViRest_fs ...
                                  + HL.ge .* Psi_ei_fs(Vi_fs_grid) .* Phi_ei ...      %E-to-I
                                  + HL.gi * Psi_ii_fs(Vi_fs_grid) .* Phi_ii_fs ...      %I-to-I
                                  + 0 .* (lap_synap * Vi_fs_grid));

        % 4. update the firing rates
        sig_offset = 18;
        Qe_grid = HL.Qe_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_e) .* (Ve_grid - HL.theta_e)))) ...     % The E voltage must be big enough,
                  - HL.Qe_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_e/5) .* (Ve_grid - (HL.theta_e + sig_offset)))));     % ... but not too big.
        Qi_grid = HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i) .* (Vi_grid - HL.theta_i)))) ...     % The I voltage must be big enough,
                  - HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i/5) .* (Vi_grid - (HL.theta_i + sig_offset)))));     % ... but not too big.
        Qi_fs_grid = HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i) .* (Vi_fs_grid - HL.theta_i)))) ...     % The I voltage must be big enough,
                   - HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i/5) .* (Vi_fs_grid - (HL.theta_i + sig_offset)))));
        
        Qe_grid = max(Qe_grid, 0);
        Qi_grid = max(Qi_grid, 0);
        Qi_fs_grid = max(Qi_fs_grid, 0);
        
        % Qe_grid(locs(:,3) > -9.5 & locs(:,3) < -9.3) = 0;
        % Qi_fs_grid(map > 0) = 0;
              
        % 5. update extracellular ion
        joint_Q = Qi_fs_grid + HL.prodRatio * (Qe_grid + Qi_grid);
        % joint_Q = Qe_grid + Qi_grid + Qi_fs_grid;
        K_1 = K + dt/HL.tau_K * (-HL.k_decay .* K ...   % decay term.
                + HL.kR .* joint_Q ./ (1+exp(-(joint_Q - 3))) ... % reaction term.
                + HL.kD * (lap_diff * K));          % diffusion term.

        % 6. update inhibitory gap junction strength, and resting voltages
        D22_1         = D22        + dt/HL.tau_dD * (HL.KtoD * (K-5));
        
        PNa = 0.04; NaO = 145; NaI = 10;
        PCl = 0.45; ClO = 110; ClI = 6;
        KI = 155; Kinit = 5;
        
        ghk = 26.7123 * log((K + PNa*NaO + PCl*ClI) ./ (KI + PNa*NaI + PCl*ClO));
        HL.Ve_rest = ghk;
        HL.Vi_rest = ghk;
        % del_VeRest_1  = ghk - HL.Ve_rest;
        % del_ViRest_1  = ghk - HL.Vi_rest;
        % del_VeRest_1  = del_VeRest + dt/HL.tau_dVe * (HL.KtoVe*K);
        % del_ViRest_1  = del_ViRest + dt/HL.tau_dVi * (HL.KtoVi*K);
        
        PK = 1.3;
        ghk = 26.7123 * log((PK*K + PNa*NaO + PCl*ClI) ./ (PK*(KI - 18*(K-Kinit)) + PNa*NaI + PCl*ClO));
        HL.Vi_rest_fs = ghk;
        % del_ViRest_fs_1  = ghk - HL.Vi_rest;
        % del_ViRest_fs_1  = 6 * K;
        % del_ViRest_fs_1  = del_ViRest_fs + dt/HL.tau_dVi * (HL.KtoVi_fs*K);

        % 7. update dynamic variables
        phi2_ee = phi2_ee_1;
        phi_ee  = phi_ee_1;
        phi2_ei = phi2_ei_1;
        phi_ei  = phi_ei_1;

        F_ee   = F_ee_1;
        Phi_ee = Phi_ee_1;
        F_ei   = F_ei_1;
        Phi_ei = Phi_ei_1;

        F_ie   = F_ie_1;
        Phi_ie = Phi_ie_1;
        F_ii   = F_ii_1;
        F_ii_fs   = F_ii_fs_1;
        Phi_ii = Phi_ii_1;
        Phi_ii_fs = Phi_ii_fs_1;

        Ve_grid = Ve_grid_1;
        Vi_grid = Vi_grid_1;
        Vi_fs_grid = Vi_fs_grid_1;
        D22 = max(D22_1, HL.D22min);                  % the inhibitory gap junctions cannot pass below a minimum value of 0.1.
        D11 = D22 / 100;                              % see definition in [Steyn-Ross et al PRX 2013, Table I].
        
        % del_VeRest = min(del_VeRest_1, 30);            % the excitatory population resting voltage cannot pass above a maximum value of 1.5.    
        if ~isnan(source_del_VeRest)
            % Qe_grid(map > 0) = source_del_VeRest;
            % del_VeRest(map > 0) = source_del_VeRest * map(map > 0); % + del_VeRest(map > 0);     % set the "source" locations' excitatory population resting voltage
        end
  
        % del_ViRest = min(del_ViRest_1, 30);           % the inhibitory population resting voltage cannot pass above a maximum value of 0.8.
        % del_ViRest_fs = min(del_ViRest_fs_1, 60);
        K = max(min(K_1, 12), 5);                               % the extracellular ion cannot pass above a maximum value of 1.0.

        % sanity check!
        if any(any(isnan(Qe_grid)))
          error('Sigmoid generated NaNs!! (Either increase dx or reduce dt)');
        end
    end

    % save the values at the end of time period
    last.Qe = Qe_grid;
    last.Qi = Qi_grid;
    last.Qi_fs = Qi_fs_grid;
    last.Ve = Ve_grid;
    last.Vi = Vi_grid;
    last.Vi_fs = Vi_fs_grid;

    last.phi_ee = phi_ee;
    last.phi_ei = phi_ei;
    last.phi2_ee = phi2_ee;
    last.phi2_ei = phi2_ei;

    last.F_ee = F_ee;
    last.F_ei = F_ei;
    last.F_ie = F_ie;
    last.F_ii = F_ii;
    last.F_ii_fs = F_ii_fs;

    last.Phi_ee = Phi_ee;
    last.Phi_ei = Phi_ei;
    last.Phi_ie = Phi_ie;
    last.Phi_ii = Phi_ii;
    last.Phi_ii_fs = Phi_ii_fs;

    last.D11 = D11;
    last.D22 = D22;
    last.dVe = del_VeRest;
    last.dVi = del_ViRest;
    last.dVi_fs = del_ViRest_fs;
    last.K = K;
end

%------------------------------------------------------------------------
function weight = Psi_ee(V)
    % e-to-e reversal-potential weighting function

    global HL
    weight = (HL.Ve_rev - V)./(HL.Ve_rev - HL.Ve_rest);
end

%------------------------------------------------------------------------
function weight = Psi_ei(V)
    % e-to-i reversal-potential weighting function

    global HL
    weight = (HL.Ve_rev - V)./(HL.Ve_rev - HL.Vi_rest);
end
function weight = Psi_ei_fs(V)
    % e-to-i reversal-potential weighting function

    global HL
    weight = (HL.Ve_rev - V)./(HL.Ve_rev - HL.Vi_rest_fs);
end

%------------------------------------------------------------------------
function weight = Psi_ie(V)
    % i-to-e reversal-potential weighting function

    global HL
    weight = (HL.Vi_rev - V)./(HL.Vi_rev - HL.Ve_rest);
end

%------------------------------------------------------------------------
function weight = Psi_ii(V)
    % i-to-i reversal potential weighting function

    global HL
    weight = (HL.Vi_rev - V)./(HL.Vi_rev - HL.Vi_rest);
end
function weight = Psi_ii_fs(V)
    % i-to-i reversal potential weighting function

    global HL
    weight = (HL.Vi_rev - V)./(HL.Vi_rev - HL.Vi_rest_fs);
end

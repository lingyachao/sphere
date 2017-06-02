function [time,last,Qe_micro,Qe_macro] = seizing_cortical_field_edge(source_del_VeRest, map, time_end, IC, ...
    laplacian, avg_D, edge, avg_flux, micro_idx, macro_idx)
 
    % noise level
    noise = 0.5;             

    % parameters for proportion of extracellular potassium.
    tau_K = 200;    % time-constant (/s).
    k_decay = 0.1;  % decay rate (/s).
    kD = 1;         % diffusion coefficient (cm^2/s).
    KtoVe = 10;     % impact on excitatory population resting voltage.
    KtoVi = 10;     % impact on inhibitory population resting voltage.
    KtoD  = -50;    % impact on inhibitory gap junction strength.
    kR    = 0.15;   % scale reaction term. 

    tau_dD  = 200;  % inhibitory gap junction time-constant (/s).
    tau_dVe = 250;  % excitatory population resting voltage time-constant (/s).
    tau_dVi = 250;  % inhibitory population resting voltage time-constant (/s).

    % set no. of sampling points (must be even!) along each axis of cortical grid
    N = size(laplacian, 1);

    % initialize constants for Seizing Cortical Model
    global HL
    HL = SCM_init_globs;

    % define synaptic strengths
    rho_e = HL.ge;
    rho_i = HL.gi;

    % initialize random number generator (from input argument)
    rand_state = sum(100*clock);
    rng(rand_state); % randn('state', rand_state);

    noise_sf = 0.2*20*noise;    % noise scale-factor
    noise_sc = 0.2;             % subcortical noise

    HL.v = 280;                 % axonal conduction velocity (cm/s), [original = 140 cm/s]
    HL.Lambda = 4.0;			% inverse-length scale for connectivity (/cm)
    HL.gamma_e = 170;           % EPSP decay rate (/s)
    HL.gamma_i = 50;            % IPSP decay rate (/s)
    HL.tau_e = 0.02;			% excit neuron time-constant (/s) [original = 0.04 s]
    HL.tau_i = 0.02;			% inhib neuron time-constant (/s) [original = 0.04 s]

    % set time resolution
    dt = 0.2 * 1e-3;

    % number of time-steps for simulation
    Nsteps = round(time_end/dt);
    time   = (0:Nsteps-1)' * dt;

    % use as initial conditions the "last" values of previous simulation.
    Qe_grid = IC.Qe;
    Qi_grid = IC.Qi;
    Ve_grid = IC.Ve;
    Vi_grid = IC.Vi;

    phi_ee = IC.phi_ee;
    phi_ei = IC.phi_ei;
    phi2_ee = IC.phi2_ee;
    phi2_ei = IC.phi2_ei;

    F_ee = IC.F_ee;
    F_ei = IC.F_ei;
    F_ie = IC.F_ie;
    F_ii = IC.F_ii;

    Phi_ee = IC.Phi_ee;
    Phi_ei = IC.Phi_ei;
    Phi_ie = IC.Phi_ie;
    Phi_ii = IC.Phi_ii;

    D11 = IC.D11;
    D22 = IC.D22;
    del_VeRest = IC.dVe;
    del_ViRest = IC.dVi;
    K = IC.K;

    % allocate matrices for recording Qe (fine time scale) at micro/macro
    % nodes
    Qe_micro = zeros(Nsteps, length(micro_idx));
    Qe_macro = zeros(Nsteps, length(macro_idx));
    
    % noise-amplitude coefficients for subcortical flux (note 1/sqrt(dt) factor)
    B_ee = noise_sf * sqrt(noise_sc * HL.phi_ee_sc / dt);
    B_ei = noise_sf * sqrt(noise_sc * HL.phi_ei_sc / dt);

    for i = 1: Nsteps

        % 0. record Qe
        Qe_micro(i,:) = Qe_grid(micro_idx);
        Qe_macro(i,:) = Qe_grid(macro_idx);
        
        % 1. update wave equations
        phi2_ee_1 = phi2_ee + dt * (-2 * HL.v * HL.Lambda * phi2_ee ...
                    - (HL.v * HL.Lambda)^2 * phi_ee ...
                    + (HL.v * HL.Lambda)^2 * Qe_grid)...
                    + dt * (HL.v / avg_D)^2 * (laplacian * phi_ee);
        phi_ee_1 = phi_ee + dt * phi2_ee;

        phi2_ei_1 = phi2_ei + dt * (-2 * HL.v * HL.Lambda * phi2_ei ...
                    - (HL.v * HL.Lambda)^2 * phi_ei ...
                    + (HL.v * HL.Lambda)^2 * Qe_grid) ...
                    + dt * (HL.v / avg_D)^2 * (laplacian * phi_ei);
        phi_ei_1 = phi_ei + dt*phi2_ei;

        % 2. update the 4 synaptic flux equations (include sc noise)

        %%%% E-to-E %%%%
        F_ee_1   = F_ee + dt * HL.gamma_e^2 * (-2/HL.gamma_e*F_ee - Phi_ee ...
                        + HL.Nee_a * phi_ee ...         % long range
                        + HL.Nee_b * Qe_grid ...        % short range
                        + noise_sc * HL.phi_ee_sc ...   % subcortical (tonic)
                        + B_ee * randn(N, 1));          % subcortical (random)
        Phi_ee_1 = Phi_ee + dt*F_ee;

        %%%% E-to-I %%%%
        F_ei_1   = F_ei + dt * HL.gamma_e^2 * (-2/HL.gamma_e*F_ei - Phi_ei ...
                        + HL.Nei_a * phi_ei ...         % long range
                        + HL.Nei_b * Qe_grid ...        % short range
                        + noise_sc * HL.phi_ei_sc ...   % subcortical (tonic)
                        + B_ei * randn(N, 1));          % subcortical (random)
        Phi_ei_1 = Phi_ei + dt*F_ei;

        %%%% I-to-E %%%%
        F_ie_1   = F_ie + dt * HL.gamma_i^2 * (-2/HL.gamma_i*F_ie - Phi_ie ...
                        + HL.Nie_b * Qi_grid);          % short range
        Phi_ie_1 = Phi_ie + dt*F_ie;

        %%%% I-to-I %%%%
        F_ii_1   = F_ii + dt * HL.gamma_i^2 * (-2/HL.gamma_i*F_ii - Phi_ii ...
                        + HL.Nii_b * Qi_grid);          % short range
        Phi_ii_1 = Phi_ii + dt*F_ii;

        % 3. update the soma voltages

        Ve_grid_1 = Ve_grid + dt/HL.tau_e * ((HL.Ve_rest - Ve_grid) + del_VeRest ...
              + rho_e * Psi_ee(Ve_grid) .* Phi_ee ...      %E-to-E
              + rho_i * Psi_ie(Ve_grid) .* Phi_ie ...      %I-to-E
              + D11 .* (laplacian * Ve_grid));

        Vi_grid_1 = Vi_grid + dt/HL.tau_i * ((HL.Vi_rest - Vi_grid) + del_ViRest ...
              + rho_e * Psi_ei(Vi_grid) .* Phi_ei ...      %E-to-I
              + rho_i * Psi_ii(Vi_grid) .* Phi_ii ...      %I-to-I
              + D22 .* (laplacian * Vi_grid));

        % 4. update the firing rates
        Qe_grid = HL.Qe_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_e) .* (Ve_grid - HL.theta_e)))) ...     %The E voltage must be big enough,
                - HL.Qe_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_e) .* (Ve_grid - (HL.theta_e+30)))));   %... but not too big.
        Qi_grid = HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i) .* (Vi_grid - HL.theta_i)))) ...     %The I voltage must be big enough,
                - HL.Qi_max * (1./(1+exp(-pi/(sqrt(3)*HL.sigma_i) .* (Vi_grid - (HL.theta_i+30)))));   %... but not too big.

        % 5. update extracellular ion
        K_1 = K + dt/tau_K * (-k_decay*K ...   % decay term.                  
                + kR * (Qe_grid + Qi_grid)./(1+exp(-((Qe_grid + Qi_grid)-15))) ... % reaction term.
                + kD * (laplacian * K));     % diffusion term.

        % 6. update inhibitory gap junction strength, and resting voltages
        D22_1         = D22        + dt/tau_dD  * (KtoD*K);
        del_VeRest_1  = del_VeRest + dt/tau_dVe * (KtoVe*K);
        del_ViRest_1  = del_ViRest + dt/tau_dVi * (KtoVi*K);

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
        Phi_ii = Phi_ii_1;

        Ve_grid = Ve_grid_1;
        Vi_grid = Vi_grid_1;
        D22 = max(D22_1,0.1);                       % the inhibitory gap junctions cannot pass below a minimum value of 0.1.
        D11 = D22/100;                              % see definition in [Steyn-Ross et al PRX 2013, Table I].

        del_VeRest = min(del_VeRest_1, 1.5);        % the excitatory population resting voltage cannot pass above a maximum value of 1.5.    
        del_VeRest(map == 1) = source_del_VeRest;   % set the "source" locations' excitatory population resting voltage
        del_ViRest = min(del_ViRest_1,0.8);         % the inhibitory population resting voltage cannot pass above a maximum value of 0.8.
        K = min(K_1,1);                             % the extracellular ion cannot pass above a maximum value of 1.0.

        % 8. no flux boundary conditions
        phi2_ee(edge) = avg_flux * phi2_ee;
        phi_ee(edge)  = avg_flux * phi_ee;
        phi2_ei(edge) = avg_flux * phi2_ei;
        phi_ei(edge)  = avg_flux * phi_ei;

        F_ee(edge)   = avg_flux * F_ee;
        Phi_ee(edge) = avg_flux * Phi_ee;
        F_ei(edge)   = avg_flux * F_ei;
        Phi_ei(edge) = avg_flux * Phi_ei;

        F_ie(edge)   = avg_flux * F_ie;
        Phi_ie(edge) = avg_flux * Phi_ie;
        F_ii(edge)   = avg_flux * F_ii;
        Phi_ii(edge) = avg_flux * Phi_ii;

        Ve_grid(edge) = avg_flux * Ve_grid;
        Vi_grid(edge) = avg_flux * Vi_grid;
        D22(edge) = avg_flux * D22;
        D11(edge) = avg_flux * D11;

        del_VeRest(edge) = avg_flux * del_VeRest;    
        del_ViRest(edge) = avg_flux * del_ViRest;
        K(edge) = avg_flux * K;
        Qe_grid(edge) = avg_flux * Qe_grid;
        Qi_grid(edge) = avg_flux * Qi_grid;
        
        % sanity check!
        if any(any(isnan(Qe_grid)))
          error('Sigmoid generated NaNs!! (Either increase dx or reduce dt)');
        end
    end

    % save the values at the end of time period
    last.Qe = Qe_grid;
    last.Qi = Qi_grid;
    last.Ve = Ve_grid;
    last.Vi = Vi_grid;

    last.phi_ee = phi_ee;
    last.phi_ei = phi_ei;
    last.phi2_ee = phi2_ee;
    last.phi2_ei = phi2_ei;

    last.F_ee = F_ee;
    last.F_ei = F_ei;
    last.F_ie = F_ie;
    last.F_ii = F_ii;

    last.Phi_ee = Phi_ee;
    last.Phi_ei = Phi_ei;
    last.Phi_ie = Phi_ie;
    last.Phi_ii = Phi_ii;

    last.D11 = D11;
    last.D22 = D22;
    last.dVe = del_VeRest;
    last.dVi = del_ViRest;
    last.K = K;

end

%------------------------------------------------------------------------
function weight = Psi_ee(V)
    % e-to-e reversal-potential weighting function

    global HL
    weight = (HL.Ve_rev - V)/(HL.Ve_rev - HL.Ve_rest);
end

%------------------------------------------------------------------------
function weight = Psi_ei(V)
    % e-to-i reversal-potential weighting function

    global HL
    weight = (HL.Ve_rev - V)/(HL.Ve_rev - HL.Vi_rest);
end

%------------------------------------------------------------------------
function weight = Psi_ie(V)
    % i-to-e reversal-potential weighting function

    global HL
    weight = (HL.Vi_rev - V)/(HL.Vi_rev - HL.Ve_rest);
end

%------------------------------------------------------------------------
function weight = Psi_ii(V)
    % i-to-i reversal potential weighting function

    global HL
    weight = (HL.Vi_rev - V)/(HL.Vi_rev - HL.Vi_rest);
end

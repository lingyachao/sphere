function H = SCM_init_globs(N)

    % time constants
    H.v = 280;                                      % axonal conduction velocity (cm/s), [original = 140 cm/s]
    H.Lambda = 4.0;                                 % inverse-length scale for connectivity (/cm)
    H.gamma_e = 170;                                % EPSP decay rate (/s)
    H.gamma_i = 50;                                 % IPSP decay rate (/s)
    H.tau_e = 0.02;                                 % excit neuron time-constant (/s) [original = 0.04 s]
    H.tau_i = 0.02;                                 % inhib neuron time-constant (/s) [original = 0.04 s]

    % noise factors
    H.noise = 0.5;
    % H.noise = 1.7; % to generate spontaneous seizures, use this and set
    % Qi-fs at focus to 0 and last.dVe(:) = 0.6;
    H.noise_sf = 0.2*20*H.noise;                    % noise scale-factor
    H.noise_sc = 0.2;                               % subcortical noise
    
    % sigmoid characteristics
    [H.Qe_max,  H.Qi_max]  = deal(30,    60);       % sigmoid maximum (s^-1)
    [H.theta_e, H.theta_i] = deal(-58.5, -58.5);	% sigmoid threshold (mV)
    [H.sigma_e, H.sigma_i] = deal(3.0,   5.0);		% sigmoid 'width' (mV)

    % gain per synapse at resting voltage (millivolt.sec)
    [H.ge, H.gi] = deal(1.00e-3, -1.05e-3);

    % voltage limits
    [H.Ve_rev,  H.Vi_rev]  = deal(0, -70);          % reversal potential (mV)
    [H.Ve_rest, H.Vi_rest] = deal(-64, -64);        % resting potential (mV)

    % connectivities: j-->k convention (dimensionless)
    [H.Nee_a, H.Nei_a] = deal(2000, 2000);			% cortico-cortical
    [H.Nee_b, H.Nei_b] = deal(800, 800);
    [H.Nie_b, H.Nii_b] = deal(600, 600);
    [H.Nie_fs, H.Nii_fs] = deal(300, 100);
    [H.Nee_sc,H.Nei_sc]= deal(50, 50);              % subcortical

    % parameters for proportion of extracellular potassium.
    H.tau_K = 200;                                  % time-constant (/s).
    H.k_decay = 0.1;                                % decay rate (/s).
    H.kR = 0.15;                                    % scale reaction term. 
    H.kD = 1;                                       % diffusion coefficient (cm^2/s).
    
    H.KtoVe = 10;                                   % impact on excitatory population resting voltage.
    H.KtoVi = 10;                                   % impact on inhibitory population resting voltage.
    H.KtoD  = -50;                                  % impact on inhibitory gap junction strength.

    H.tau_dVe = 250;                                % excitatory population resting voltage time-constant (/s).
    H.tau_dVi = 250;                                % inhibitory population resting voltage time-constant (/s).
    H.tau_dD  = 200;                                % inhibitory gap junction time-constant (/s).
    
    % Nee and Nie totals for cortico-cortical plus intracortical
    H.Nee_ab = H.Nee_a + H.Nee_b;
    H.Nei_ab = H.Nei_a + H.Nei_b;

    % default subcortical fluxes
    H.phi_ee_sc = H.Nee_sc * H.Qe_max;
    H.phi_ei_sc = H.Nei_sc * H.Qe_max;

    % d/dV derivatives of psi_ij weighting functions
    H.d_psi_ee = -1/(H.Ve_rev - H.Ve_rest);
    H.d_psi_ei = -1/(H.Ve_rev - H.Vi_rest);
    H.d_psi_ie = -1/(H.Vi_rev - H.Ve_rest);
    H.d_psi_ii = -1/(H.Vi_rev - H.Vi_rest);

    % non-homogenous parameters
    H.Nie_b = H.Nie_b * ones(N, 1);
    H.Nii_b = H.Nii_b * ones(N, 1);
    H.Vi_rest = H.Vi_rest * ones(N, 1);
    H.ge = H.ge * ones(N, 1);
    H.phi_ee_sc = H.phi_ee_sc * ones(N, 1);
    
    return
end
%------------------------------------------------------------------------

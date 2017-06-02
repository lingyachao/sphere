function IC = make_IC(N)
% output an initial state

    load('seizing_cortical_field_IC.mat');

    IC.Qe = fill(N, last.Qe);
    IC.Qi = fill(N, last.Qi);
    IC.Ve = fill(N, last.Ve);
    IC.Vi = fill(N, last.Vi);

    IC.phi_ee = fill(N, last.phi_ee);
    IC.phi_ei = fill(N, last.phi_ei);
    IC.phi2_ee = fill(N, last.phi2_ee);
    IC.phi2_ei = fill(N, last.phi2_ei);

    IC.F_ee = fill(N, last.F_ee);
    IC.F_ei = fill(N, last.F_ei);
    IC.F_ie = fill(N, last.F_ie);
    IC.F_ii = fill(N, last.F_ii);

    IC.Phi_ee = fill(N, last.Phi_ee);
    IC.Phi_ei = fill(N, last.Phi_ei);
    IC.Phi_ie = fill(N, last.Phi_ie);
    IC.Phi_ii = fill(N, last.Phi_ii);

    IC.D11 = fill(N, last.D11);
    IC.D22 = fill(N, last.D22);
    IC.dVe = fill(N, last.dVe);
    IC.dVi = fill(N, last.dVi);
    IC.K = fill(N, last.K);
    
end

function new = fill(N, old)

    N0 = length(old(:));
    q = floor(N/N0);
    r = rem(N, N0);
    new = [repmat(old(:), q, 1); old(1:r)'];

end



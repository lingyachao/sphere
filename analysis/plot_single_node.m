


for k = 1:K
    fprintf(['Read in ' num2str(k) '\n']);
    load([RAW_DIR 'seizing_cortical_field_k_'  num2str(k) '.mat'], 'last');
    
    Qe_1(k) = last.Qe(node_id);
    Qi_1(k) = last.Qi(node_id);
    Ve_1(k) = last.Ve(node_id);
    Vi_1(k) = last.Vi(node_id);
    D22_1(k) = last.D22(node_id);
    dVe_1(k) = last.dVe(node_id);
    dVi_1(k) = last.dVi(node_id);
    K_1(k) = last.K(node_id);
end

single_node = table(Qe_1, Qi_1, Ve_1, Vi_1, D22_1, dVe_1, dVi_1, K_1);
plot((1:K) * T0, table2array(single_node));
legend('Qe', 'Qi', 'Ve', 'Vi', 'D22', 'dVe', 'dVi', 'K');
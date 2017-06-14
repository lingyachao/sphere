for k=1:200
   
    load(['./data/brain_N40962_06131715_full_data/raw/seizing_cortical_field_k_' num2str(k) '.mat']);
    
    fine = rmfield(fine, 'Qe_micro');
    fine = rmfield(fine, 'Ve_micro');
    fine = rmfield(fine, 'Qe_focus');
    fine = rmfield(fine, 'Ve_focus');
    fine.Qe_lessihb = fine.Qe_macro;
    fine.Ve_lessihb = fine.Ve_macro;
    fine = rmfield(fine, 'Qe_macro');
    fine = rmfield(fine, 'Ve_macro');
    
    save(['./data/brain_N40962_06141205_full_data/raw/seizing_cortical_field_k_' num2str(k) '.mat'], ...
            'samp_time', 'last', 'fine');
end    
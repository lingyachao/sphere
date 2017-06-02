global HL
HL = SCM_init_globs;

HL.kR = 0;
HL.Nie_b = HL.Nie_b * ones(N, 1);
HL.Nii_b = HL.Nii_b * ones(N, 1);
HL.Vi_rest = HL.Vi_rest * ones(N, 1);

ge_ticks = 0.3:0.02:1;
phi_ticks = 12:45;

V_diff = NaN(length(phi_ticks), length(ge_ticks));

for i = 1:length(ge_ticks)
    for j = 1:length(phi_ticks)
        HL.ge = 1.00e-3 * ge_ticks(i);
        HL.phi_ee_sc = 1500 * phi_ticks(j);
        run_seizing_cortical_field;
        end_samp = Ve_samp(end*0.8:end);
        V_diff(j,i) = max(end_samp) - min(end_samp);
    end
end

imagesc(V_diff);
set(gca,'YDir','normal');
set(gca,'XTick', 1:length(ge_ticks));
set(gca,'YTick', 1:length(phi_ticks));
set(gca,'XTickLabel', cellstr(num2str(ge_ticks')));
set(gca,'YTickLabel', cellstr(num2str(phi_ticks')));
colorbar;


N = 10242;

global HL
HL = SCM_init_globs;

HL.kR = 0;

% ge_ticks = 0.3:0.02:1;
% phi_ticks = 12:45;
ge_ticks = 0.8:0.1:1.5;
phi_ticks = 1:5;

V_diff = NaN(length(phi_ticks), length(ge_ticks));

focus_indices = 1:7;

for i = 1:length(ge_ticks)
    for j = 1:length(phi_ticks)
        
        fprintf(['Running ' num2str(ge_ticks(i)) ' ' num2str(phi_ticks(j)) '\n']);
        
        HL.ge(focus_indices) = 1.00e-3 * ge_ticks(i);
        HL.phi_ee_sc(focus_indices) = 1500 * phi_ticks(j);
        run_seizing_cortical_field;
        % end_samp = Ve_samp(end*0.8:end);
        % V_diff(j,i) = max(end_samp) - min(end_samp);
        V_diff(j,i) = max(Ve_samp);
    end
end

imagesc(V_diff);
set(gca,'YDir','normal');
set(gca,'XTick', 1:length(ge_ticks));
set(gca,'YTick', 1:length(phi_ticks));
set(gca,'XTickLabel', cellstr(num2str(ge_ticks')));
set(gca,'YTickLabel', cellstr(num2str(phi_ticks')));
colorbar;


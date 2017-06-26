N = 42;

global source_drive
global coupling

Ve_samp_focus = NaN(1500, 10, 5);
Ve_samp_homo = NaN(1500, 10, 5);

for i = 1:10
    for j = 1:5
        
        source_drive = i;
        coupling = j * 0.5;
        
        subplot(5, 10, i + 10*(j-1));
        
        lessihb_filter = locs(:,3) < -6;
        run_seizing_cortical_field;
        Ve_samp_focus(:,i,j) = Ve_samp;
        
        hold on;
        
        lessihb_filter = true(N, 1);
        run_seizing_cortical_field;
        Ve_samp_homo(:,i,j) = Ve_samp;
        
        ylim([-70 -45]);
        title(['dVe = ' num2str(source_drive) ' D22 = ' num2str(coupling)]);
    end
end
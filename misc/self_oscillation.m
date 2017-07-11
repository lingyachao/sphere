global source_drive
global coupling

Ve_samp_focus = NaN(1500, 8, 5);
Vi_samp_focus = NaN(1500, 8, 5);
Ve_samp_homo = NaN(1500, 8, 5);
Vi_samp_homo = NaN(1500, 8, 5);

for i = 1:8
    for j = 1:5

        source_drive = i;
        coupling = j * 0.5;
        fprintf(['dVe = ' num2str(source_drive) ' D22 = ' num2str(coupling) '\n']);
        
        subplot(5, 8, i + 8*(j-1));
        
        N = 10242;
        run_seizing_cortical_field;
        Ve_samp_focus(:,i,j) = Ve_samp;
        Vi_samp_focus(:,i,j) = Vi_samp;
        
        hold on;
        
        N = 12;
        run_seizing_cortical_field;
        Ve_samp_homo(:,i,j) = Ve_samp;
        Vi_samp_homo(:,i,j) = Vi_samp;
        
        ylim([-70 -45]);
        title(['dVe = ' num2str(source_drive) ' D22 = ' num2str(coupling)], 'FontSize', 10);
        drawnow;
    end
end

save('./images/self_oscillation_grid_12.mat', ...
    'Ve_samp_focus', 'Vi_samp_focus', ...
    'Ve_samp_homo', 'Vi_samp_homo');

savefig('./images/self_oscillation_grid_12.fig');
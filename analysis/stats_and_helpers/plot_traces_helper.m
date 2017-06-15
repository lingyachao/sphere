function plot_traces_helper(t, y, plot_title, plot_ylabel, plot_ylim)
    plot(t, y);
    % line([30 30], [-100 100], 'Color', 'r', 'LineWidth', 2);
    % line([150 150], [-100 100], 'Color', 'r', 'LineWidth', 2);
    title(plot_title);
    xlabel('time (s)');
    ylabel(plot_ylabel);
    ylim(plot_ylim);
    
    if size(y, 2) == 3
        legend('focus', 'less-inhibition', 'normal');
    else
        legend('focus', 'less-inhibition 1', 'less-inhibition 2', 'normal');
    end
end


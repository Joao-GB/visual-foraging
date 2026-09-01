function p = plotPSAPermTest(T, permT, params)
    figure('Name', 'PSA permutation test', 'Color', 'w', 'Position', [100 100 1400 600]);
    hold on;
    histogram(permT, 'Normalization', 'count', 'FaceColor', params.darkBlue, 'EdgeColor', 'w', 'FaceAlpha', 0.75);
    set(gca,'TickDir', 'out', 'Box', 'off','FontSize', 12);
    % P-valor two sided
    p = (1 + sum(abs(permT) >= abs(T))) / (numel(permT) + 1);
    xline(T, '--k', 'LineWidth', 2);
    yl = ylim;
    text(T, yl(2)*0.25, sprintf('  p = %.4g', p), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 10);
    xlabel('\Deltad'' médio intra-sujeitos');
    ylabel('Proporção');
end
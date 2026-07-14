function renderPSAforagingPerformance(plotData, trialCounts, ylabelStr, figName, plotTitle, drP, isSens)
    figure('Name', figName, 'Color', 'w', 'Position', [100 200 900 450]);
    
    t = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'normal');
    xlabelPSA = ["Sacádico", "Não-sacádico"];
    subTitles = ["Forrageamento Correto", "Forrageamento Incorreto"];
    
    for i = 1:2
        ax = nexttile;
        currentData = plotData{i};
        currentCounts = trialCounts{i};
        
        % Configura espessura com base na métrica
        if isSens
            b = bar(currentData, 0.3, 'FaceColor', 'flat', 'LineWidth', 1.5);
        else
            b = bar(currentData, 'FaceColor', 'flat');
        end
        
        b.EdgeColor = [0 0 0];
        b.CData(1,:) = drP.darkBlue;
        b.CData(2,:) = drP.paleBrown;
        
        set(ax, 'TickDir', 'out', 'Box', 'off');
        xticklabels(xlabelPSA);
        ylabel(ylabelStr); 
        title(subTitles{i});
        grid on;
        
        if isSens
            localMax = max(currentData);
            if localMax > 0
                ylim([0 localMax * 1.2]); 
            else
                ylim([0 1]);
            end
        else
            ylim([0 100]);
        end
        text(ax, .95, .95, sprintf('n = %d', currentCounts(1)), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontSize', 9,'FontWeight', 'bold', 'Color', .4*[1 1 1]);
    end
    
    title(t, plotTitle, 'FontSize', 13, 'FontWeight', 'bold');
end
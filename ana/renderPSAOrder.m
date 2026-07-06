function renderPSAOrder(plotDataCell, xlabelOrder, ylabelOrder, figName, plotTitle, drP, isSens)
    figure('Name', figName, 'Color', 'w', 'Position', [100 100 1400 600]);
    
    % Gera um plot com 13 de largura, pois a 10a será para um risquinho
    t = tiledlayout(2, 13, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    colStarts = [1, 4, 7, 11]; 
    
    % A cada 3 há um subplot
    for r = 1:2
        for c = 1:4 
            tileIdx = (r - 1) * 13 + colStarts(c);
            ax = nexttile(tileIdx, [1, 3]); 
            if r == 1 && c == 3, ax3 = ax; end
            if r == 1 && c == 4, ax4 = ax; end

            if isSens
                % Your specific thin custom design configuration
                b = bar(plotDataCell{r, c}, 0.3, 'FaceColor', 'flat', 'LineWidth', 1.5);
                grid on; 
            else
                % Original accuracy bar thickness configurations
                b = bar(plotDataCell{r, c}, 'FaceColor', 'flat');
                grid on;
            end
            
            b.EdgeColor = [0 0 0]; 
            
            currentLabels = xlabelOrder{r, c};
            for barIdx = 1:length(currentLabels)
                switch currentLabels(barIdx)
                    case "F"
                        b.CData(barIdx, :) = [1 1 1];
                    case "S"
                        b.CData(barIdx, :) = drP.darkBlue;
                    case "N"
                        b.CData(barIdx, :) = drP.paleBrown;
                end
            end
            
            set(ax, 'TickDir', 'out', 'Box', 'off');
            xticklabels(currentLabels);

            if isSens
                % Find the maximum value present ONLY inside this specific plot's bars
                localMax = max(plotDataCell{r, c});
                if localMax > 0
                    ylim([0 localMax * 1.2]);
                else
                    ylim([0 1]); % Fallback if all data values are 0
                end
            else
                ylim([0 100]); % Keeps accuracy locked between 0% and 100%
            end
            
            if c == 1
                ylabel(ylabelOrder{r});
            end
            
            if c == 4
                set(ax, 'Color', [0.95 0.95 0.95]); 
            end
        end
    end
    
    title(t, plotTitle, 'FontSize', 14, 'FontWeight', 'bold');
    
    drawnow; 
    
    lineX = ax3.Position(1) + ax3.Position(3) + (ax4.Position(1) - (ax3.Position(1) + ax3.Position(3)))/2;
    
    tPos = t.Position;
    lineY_start = tPos(2);
    lineY_end   = tPos(2) + tPos(4);
    
    annotation('line', [lineX lineX], [lineY_start lineY_end], ...
        'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.2);
end
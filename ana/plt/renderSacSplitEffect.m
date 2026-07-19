function renderSacSplitEffect(plotData, countData, colTitles, ylabelStr, figName, plotTitle, drP, isSens)
    nCols = length(colTitles);
    
    % A figura ajusta a largura dinamicamente baseada no número de colunas
    figWidth = max(600, nCols * 140);
    figure('Name', figName, 'Color', 'w', 'Position', [50 50 figWidth 500]);
    t = tiledlayout(2, nCols, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    % Busca o máximo global para unificar o eixo Y
    allVals = cellfun(@(x) max(x(:)), plotData);
    globalMax = max(allVals(:));
    if isnan(globalMax) || globalMax == 0, globalMax = 1; end
    
    rowLabels = {'Dentro', 'Fora'};
    
    for r = 1:2
        for c = 1:nCols
            ax = nexttile;
            
            data = plotData{r,c};
            count = countData{r,c};
            
            if isSens
                b = bar(1:2, data, 0.6, 'FaceColor', 'flat', 'LineWidth', 1.5);
                ylim([0 globalMax * 1.35]); % Respiro para o texto "n="
            else
                b = bar(1:2, data, 0.6, 'FaceColor', 'flat');
                ylim([0 115]); 
            end
            
            b.EdgeColor = [0 0 0];
            b.CData(1,:) = drP.darkBlue;
            b.CData(2,:) = drP.paleBrown;
            
            set(ax, 'TickDir', 'out', 'Box', 'off');
            xticks([1 2]);
            xticklabels({'S', 'N'});
            grid on;
            
            % Título da coluna (agnóstico: pode ser ângulo, altura, etc.)
            if r == 1
                title(colTitles{c}, 'FontSize', 12, 'FontWeight', 'bold');
            end
            
            % Rótulo do Eixo Y apenas na primeira coluna
            if c == 1
                ylabel({['\bf{', rowLabels{r}, '}\rm'], ylabelStr}, 'Interpreter', 'tex');
            end
            
            % Texto "n=..." sobre as barras
            if count > 0
                maxVal = max(data);
                if isnan(maxVal), maxVal = 0; end
                
                text(1.5, maxVal, sprintf('n=%d', count), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 8, ...
                    'FontWeight', 'bold', ...
                    'Color', [0.35 0.35 0.35]);
            end
        end
    end
    
    title(t, plotTitle, 'FontSize', 14, 'FontWeight', 'bold');
end
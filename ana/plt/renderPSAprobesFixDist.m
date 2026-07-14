function renderPSAprobesFixDist(plotData, countData, edges, ylabelStr, figName, plotTitle, drP, isSens)
    figure('Name', figName, 'Color', 'w', 'Position', [50 50 1000 850]);
    
    nBins = length(edges) - 1;
    gridSize = nBins + 1; % 4
    
    % Cria rótulos baseados nas bordas informadas
    labels = cell(1, gridSize);
    for i = 1:nBins
        if i == nBins && isinf(edges(i+1))
            labels{i} = sprintf('\\geq %g^o', edges(i));
        else
            labels{i} = sprintf('[%g, %g)^o', edges(i), edges(i+1));
        end
    end
    labels{gridSize} = 'Marginal';
    
    plotWidth = 3;  % tamanho de cada subplot bloco
    gapWidth  = 1;  % tamanho do espaço para a linha tracejada
    
    % A grade tem tamanho para subplots + marginais e para as divisórias
    gridDim = ((nBins + 1) * plotWidth) + gapWidth;
    t = tiledlayout(gridDim, gridDim, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    starts = [1:plotWidth:(nBins*plotWidth), (nBins*plotWidth) + gapWidth + 1];
    
    ax_core_right = []; ax_marg_right = []; ax_core_bottom = []; ax_marg_bottom = [];

    allVals = cellfun(@(x) max(x(:)), plotData);
    globalMax = max(allVals(:));
    if isnan(globalMax) || globalMax == 0, globalMax = 1; end
    
    for r = 1:gridSize
        for c = 1:gridSize
            % Se preferir deixar o canto inferior direito vazio como na outra função, 
            % basta descomentar a linha abaixo:
            if r == gridSize && c == gridSize, continue; end
            
            tileIdx = (starts(r) - 1) * gridDim + starts(c);
            ax = nexttile(tileIdx, [3, 3]); 
            
            % Captura as referências dos eixos vizinhos ao gap para desenhar as linhas
            if r == 1 && c == nBins, ax_core_right = ax;  end
            if r == 1 && c == gridSize, ax_marg_right = ax;  end
            if r == nBins && c == 1, ax_core_bottom = ax; end
            if r == gridSize && c == 1, ax_marg_bottom = ax; end
            
            data = plotData{r,c};
            count = countData{r,c};
            
            % Ajusta espessura baseado no tipo de plot
            if isSens
                b = bar(1:2, data, 0.3, 'FaceColor', 'flat', 'LineWidth', 1.5);
            else
                b = bar(1:2, data, 'FaceColor', 'flat');
            end
            
            b.EdgeColor = [0 0 0];
            b.CData(1,:) = drP.darkBlue;
            b.CData(2,:) = drP.paleBrown;
            
            set(ax, 'TickDir', 'out', 'Box', 'off');
            xticks([1 2]);
            xticklabels({'S', 'N'});
            grid on;
            
            % Fundo cinza claro para os plots marginais
            if r == gridSize || c == gridSize
                set(ax, 'Color', [0.93 0.93 0.95]);
            end
            
            % Escalonamento dinâmico do eixo Y
            if isSens
                ylim([0 globalMax * 1.2]);
            else
                ylim([0 115]);
            end
            
            % Rótulos dos eixos principais (apenas nas extremidades esquerdas/superiores)
            if c == 1
                ylabel({['\bf{Dist Sacc: ', labels{r}, '}\rm'], ylabelStr}, 'Interpreter', 'tex'); 
            end
            if r == 1
                title(['Dist N-Sacc: ', labels{c}]); 
            end
            
            % Rótulo de contagem (n) centralizado no topo
            if count > 0
                maxVal = max(data);
                if isnan(maxVal), maxVal = 0; end
                
                text(1.5, maxVal, sprintf('n=%d', count), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 9, ...
                    'FontWeight', 'bold', ...
                    'Color', [0.35 0.35 0.35]);
            end
        end
    end
    
    title(t, plotTitle, 'FontSize', 14, 'FontWeight', 'bold');
    
    drawnow; 
    tPos = t.Position;
    
    % Linha Vertical 
    lineX = ax_core_right.Position(1) + ax_core_right.Position(3) + ...
            (ax_marg_right.Position(1) - (ax_core_right.Position(1) + ax_core_right.Position(3)))/2;
            
    annotation('line', [lineX lineX], [tPos(2) tPos(2)+tPos(4)], ...
        'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.2);
        
    % Linha Horizontal
    row_core_bottom = ax_core_bottom.Position(2);
    row_marg_top    = ax_marg_bottom.Position(2) + ax_marg_bottom.Position(4);
    lineY = row_marg_top + (row_core_bottom - row_marg_top)/2;
    
    annotation('line', [tPos(1) tPos(1)+tPos(3)], [lineY lineY], ...
        'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.2);
end
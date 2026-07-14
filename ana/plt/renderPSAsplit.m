function renderPSAsplit(barMatrix, countMatrix, numBins, xlabelPSA, ylabelStr, figTitle, plotTitle, mat, isSens)
    % Ajusta proporcionalmente a largura da janela baseado no número de bins gerados
    figWidth = max(500, 250 + (numBins * 130)); 
    figure('Name', figTitle, 'Color', 'w', 'Position', [100 200 figWidth 500]);
    
    if isSens
        b = bar(barMatrix, 0.3, 'grouped', 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
    else
        b = bar(barMatrix, 'grouped', 'EdgeColor', [0 0 0]);
    end
    
    b(1).FaceColor = mat.drP.darkBlue;
    b(2).FaceColor = mat.drP.paleBrown;
    
    set(gca, 'TickDir', 'out', 'Box', 'off')
    grid on;
    
    set(gca, 'XTick', 1:numBins);
    xticklabels(xlabelPSA);
    
    ylabel(ylabelStr); 
    title(plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
    legend({'Sacádico', 'Não-sacádico'}, 'Location', 'northeast', 'Box', 'off');
    
    if isSens
        localMax = max(barMatrix(:));
        if localMax > 0
            ylim([0 localMax * 1.2]);
        else
            ylim([0 1]);
        end
    else
        ylim([0 110]);
    end

    drawnow;

    if ~isempty(countMatrix)
        for bIdx = 1:numBins
            % Encontra o valor mais alto do par para evitar que o texto cruze as barras
            maxValInCluster = max(barMatrix(bIdx, :));
            if isnan(maxValInCluster), maxValInCluster = 0; end
            
            count = countMatrix(bIdx, 1); 
            
            text(bIdx, maxValInCluster, sprintf('n = %d', count), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 9.5, ...
                'FontWeight', 'bold', ...
                'Color', [0.35 0.35 0.35]);
        end
    end
end
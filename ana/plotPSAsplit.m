function plotPSAsplit(splitValues, metricVector, trl, mat, figTitle, plotTitle)

    edges = sort(splitValues(:)'); %#ok<TRSRT> 
    numSplits = length(edges);
    numBins = numSplits + 1;

    % Linhas:  cada um dos grupos (e.g. longo e curto)
    % Colunas: 1 = Sacádico, 2 = Não-sacádico
    barMatrix = zeros(numBins, 2);
    xlabelPSA = cell(1, numBins);
    
    for bIdx = 1:numBins
        if bIdx == 1
            % Primeiro bin com < em vez de intervalo
            mask = metricVector < edges(1);
            xlabelPSA{bIdx} = sprintf('< %.3f', edges(1));
            
        elseif bIdx == numBins
            % Último bin com >= em vez de intervalo
            mask = metricVector >= edges(end);
            xlabelPSA{bIdx} = sprintf('\\geq %.3f', edges(end));
            
        else
            % Bins intermediários com intervalos
            mask = metricVector >= edges(bIdx-1) & metricVector < edges(bIdx);
            xlabelPSA{bIdx} = sprintf('[%.3f, %.3f)', edges(bIdx-1), edges(bIdx));
        end
        
        % se houver dados, calcula counts e porcentagens
        if any(mask)
            [~, countsSub] = getPSAeffect(trl(mask));
            
            % 2 para sacádico, 3 para não sacádico
            pct_s = (countsSub(1,2) / countsSub(2,2)) * 100;
            pct_n = (countsSub(1,3) / countsSub(2,3)) * 100;
            
            barMatrix(bIdx, :) = [pct_s, pct_n];
        else
            barMatrix(bIdx, :) = [0, 0];
        end
    end
    
    % Plota todos os pares de barras numa mesma figura
    figWidth = max(500, 250 + (numBins * 130)); 
    figure('Name', figTitle, 'Color', 'w', 'Position', [100 200 figWidth 500]);
    
    b = bar(barMatrix, 'grouped', 'EdgeColor', [0 0 0]);
    
    b(1).FaceColor = mat.drP.darkBlue;
    b(2).FaceColor = mat.drP.paleBrown;
    
    set(gca, 'TickDir', 'out', 'Box', 'off')
    grid on;
    
    set(gca, 'XTick', 1:numBins);
    xticklabels(xlabelPSA);
    
    ylabel('Acertos (%)'); 
    ylim([0 100]);
    title(plotTitle);
    
    legend({'Sacádico', 'Não-sacádico'}, 'Location', 'northeast', 'Box', 'off');
end
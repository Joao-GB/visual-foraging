function plotSacRangeEffect(rotProbeFix, rotProbePos, rotNSaccProbePos, sHit, nHit, drP)
    heights = [1.5, 2.5, 5, 7.5, 10, 15, 20];
    nHeights = length(heights);
    
    % Para um retângulo de largura infinita e altura H centrado no eixo X,
    % basta checar se o valor absoluto da coordenada Y do alvo é <= H/2
    absY = abs(rotNSaccProbePos(:,2));
    
    plotDataAcc = cell(2, nHeights);
    plotDataSens = cell(2, nHeights);
    countData = cell(2, nHeights);
    colTitles = cell(1, nHeights);
    
    calcDPrime = @(hit, total) norminv(min(max(hit/total, 0.001), 0.999)) * sqrt(2);
    
    for c = 1:nHeights
        colTitles{c} = sprintf('h = %g dva', heights(c)); % Formata o título
        
        halfHeight = heights(c) / 2;
        maskIn = (absY <= halfHeight);
        maskOut = ~maskIn;
        masks = {maskIn, maskOut};
        
        for r = 1:2
            m = masks{r};
            totalRows = sum(m);
            countData{r,c} = totalRows;
            
            if totalRows > 0
                plotDataAcc{r,c}  = [(sum(sHit(m)) / totalRows) * 100, (sum(nHit(m)) / totalRows) * 100];
                plotDataSens{r,c} = [calcDPrime(sum(sHit(m)), totalRows), calcDPrime(sum(nHit(m)), totalRows)];
            else
                plotDataAcc{r,c}  = [0, 0];
                plotDataSens{r,c} = [0, 0];
            end
        end
    end
    
    titleBase = 'Efeito de proximidade horizontal ao redor da direção da sacada)';
    
    renderSacSplitEffect(plotDataAcc, countData, colTitles, 'Acertos (%)', ...
        'SacRange Effect - Accuracy', [titleBase, ' (acurácia)'], drP, false);
        
    renderSacSplitEffect(plotDataSens, countData, colTitles, 'Sensibilidade (d'')', ...
        'SacRange Effect - Sensitivity', [titleBase, ' (sensibilidade)'], drP, true);
end
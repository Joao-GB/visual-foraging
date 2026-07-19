function plotSacDirEffect(rotProbeFix, rotProbePos, rotNSaccProbePos, sHit, nHit, drP)
    angles = [15, 30, 60, 90, 120, 180, 240, 270, 300, 330];
    nAngles = length(angles);
    
    theta = atan2d(rotNSaccProbePos(:,2), rotNSaccProbePos(:,1));
    absTheta = abs(theta);
    
    plotDataAcc = cell(2, nAngles);
    plotDataSens = cell(2, nAngles);
    countData = cell(2, nAngles);
    colTitles = cell(1, nAngles);
    
    calcDPrime = @(hit, total) norminv(min(max(hit/total, 0.001), 0.999)) * sqrt(2);
    
    for c = 1:nAngles
        colTitles{c} = sprintf('%g^o', angles(c)); % Formata o título
        
        halfAngle = angles(c) / 2;
        maskIn = (absTheta <= halfAngle);
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
    
    titleBase = 'Efeito da direção angular da sacada';
    
    renderSacSplitEffect(plotDataAcc, countData, colTitles, 'Acertos (%)', ...
        'SacDir Effect - Accuracy', [titleBase, ' (acurácia)'], drP, false);
        
    renderSacSplitEffect(plotDataSens, countData, colTitles, 'Sensibilidade (d'')', ...
        'SacDir Effect - Sensitivity', [titleBase, ' (sensibilidade)'], drP, true);
end
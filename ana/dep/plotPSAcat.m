function plotPSAcat(fixCat, sCat, sHit, nHit, drP)
    fixCat = fixCat(:); sCat = sCat(:); sHit = sHit(:); nHit = nHit(:);
    
    plotDataAcc = cell(3,3);
    plotDataSens = cell(3,3);
    
    % Função anônima para calcular d-prime de forma robusta (evita valores infinitos)
    calcDPrime = @(hit, total) norminv(min(max(hit/total, 0.001), 0.999)) * sqrt(2);
    
    % Calcula o plot principal 2x2
    for r = 1:2 % linhas: categoria do pré-probe (1 = alvo, 2 = distrator)
        fVal = 2 - r; % mapeia 1 -> fixCat=1, 2 -> fixCat=0
        
        for c = 1:2 % colunas: categoria do probe (1 = alvo, 2 = distrator)
            sVal = 2 - c; % 1 -> sCat=1, 2 -> sCat=0
            
            mask = (fixCat == fVal & sCat == sVal);
            
            if any(mask)
                totalRows = sum(mask);
                % Acurácia (%)
                plotDataAcc{r,c}  = [(sum(sHit(mask)) / totalRows) * 100, (sum(nHit(mask)) / totalRows) * 100];
                % Sensibilidade (d')
                plotDataSens{r,c} = [calcDPrime(sum(sHit(mask)), totalRows), calcDPrime(sum(nHit(mask)), totalRows)];
            else
                plotDataAcc{r,c}  = [0, 0];
                plotDataSens{r,c} = [0, 0];
            end
        end
    end
    
    % Calcula marginais das colunas (i.e., a última linha)
    for c = 1:2
        sVal = 2 - c;
        mask = (sCat == sVal);
        if any(mask)
            totalRows = sum(mask);
            plotDataAcc{3,c}  = [(sum(sHit(mask))/totalRows)*100, (sum(nHit(mask))/totalRows)*100];
            plotDataSens{3,c} = [calcDPrime(sum(sHit(mask)), totalRows), calcDPrime(sum(nHit(mask)), totalRows)];
        else
            plotDataAcc{3,c}  = [0, 0];
            plotDataSens{3,c} = [0, 0];
        end
    end
    
    % Calcula marginais das linhas (i.e., última coluna)
    for r = 1:2
        fVal = 2 - r;
        mask = (fixCat == fVal);
        if any(mask)
            totalRows = sum(mask);
            plotDataAcc{r,3}  = [(sum(sHit(mask))/totalRows)*100, (sum(nHit(mask))/totalRows)*100];
            plotDataSens{r,3} = [calcDPrime(sum(sHit(mask)), totalRows), calcDPrime(sum(nHit(mask)), totalRows)];
        else
            plotDataAcc{r,3}  = [0, 0];
            plotDataSens{r,3} = [0, 0];
        end
    end
    
    titleBase = 'Efeito de categoria dos estímulos';
    
    % 1. Figura de Acurácia
    renderPSAcat(plotDataAcc, 'Acertos (%)', 'PSA category effect - Accuracy', ...
                 [titleBase, ' (Acurácia)'], drP, false);
             
    % 2. Figura de Sensibilidade (d')
    renderPSAcat(plotDataSens, 'Sensibilidade (d'')', 'PSA category effect - Sensitivity', ...
                 [titleBase, ' (Sensibilidade)'], drP, true);
end
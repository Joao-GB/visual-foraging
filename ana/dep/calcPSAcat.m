function mask = calcPSAcat(fixCat, sCat, sHit, nHit)
    fixCat = fixCat(:); sCat = sCat(:); sHit = sHit(:); nHit = nHit(:);
    
    plotDataAcc = cell(3,3);
    plotDataSens = cell(3,3);
    
    % Função anônima para calcular d-prime de forma robusta (evita valores infinitos)
    calcDPrime = @(hit, total) norminv(min(max(hit/total, 0.001), 0.999)) * sqrt(2);
    mask = zeros(2,2,numel(fixCat));
    % Calcula o plot principal 2x2
    for r = 1:2 % linhas: categoria do pré-probe (1 = alvo, 2 = distrator)
        fVal = 2 - r; % mapeia 1 -> fixCat=1, 2 -> fixCat=0
        
        for c = 1:2 % colunas: categoria do probe (1 = alvo, 2 = distrator)
            sVal = 2 - c; % 1 -> sCat=1, 2 -> sCat=0
            
            mask(r,c, :) = (fixCat == fVal & sCat == sVal);
            
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
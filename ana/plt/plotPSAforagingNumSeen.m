function plotPSAforagingNumSeen(trl, drP)
    % 1. Extrair o número de itens vistos por trial
    numSeen = [trl.probeForHistIdx];
    uniqueSeen = 2:6; % Valores inteiros esperados de 2 a 6
    
    % Inicializa as matrizes para guardar os cálculos e as contagens
    barMatrixAcc  = zeros(length(uniqueSeen), 2);
    barMatrixSens = zeros(length(uniqueSeen), 2);
    countMatrix   = zeros(length(uniqueSeen), 2);
    
    % Função anônima para d-prime robusto (evita infinitos se o acerto for 0% ou 100%)
    calcDPrime = @(hit, total) norminv(min(max(hit./total, 0.001), 0.999)) * sqrt(2);
    
    % 2. Calcular o efeito PSA para cada subgrupo de numSeen
    for i = 1:length(uniqueSeen)
        currentVal = uniqueSeen(i);
        mask = (numSeen == currentVal);
        
        if any(mask)
            subTrl = trl(mask);
            [~, subCounts] = getPSAeffect(subTrl);
            
            % Registra a contagem total de trials [Sacádico, Não-sacádico]
            countMatrix(i, :) = [subCounts(2,2), subCounts(2,3)];
            
            % Acurácia (%)
            pct_s = (subCounts(1,2) / subCounts(2,2)) * 100;
            pct_n = (subCounts(1,3) / subCounts(2,3)) * 100;
            barMatrixAcc(i, :) = [pct_s, pct_n];
            
            % Sensibilidade (d')
            d_s = calcDPrime(subCounts(1,2), subCounts(2,2));
            d_n = calcDPrime(subCounts(1,3), subCounts(2,3));
            barMatrixSens(i, :) = [d_s, d_n];
        else
            barMatrixAcc(i, :)  = [0, 0];
            barMatrixSens(i, :) = [0, 0];
            countMatrix(i, :)   = [0, 0];
        end
    end
    
    titleBase = 'Efeito pré-sacádico por quantidade de estímulos vistos';
    
    % 1. Renderiza figura de Acurácia
    renderPSAforagingNumSeen(barMatrixAcc, countMatrix, uniqueSeen, 'Acertos (%)', ...
        'PSA num seen effect - Accuracy', [titleBase, ' (Acurácia)'], drP, false);
        
    % 2. Renderiza figura de Sensibilidade
    renderPSAforagingNumSeen(barMatrixSens, countMatrix, uniqueSeen, 'Sensibilidade (d'')', ...
        'PSA num seen effect - Sensitivity', [titleBase, ' (Sensibilidade)'], drP, true);
end
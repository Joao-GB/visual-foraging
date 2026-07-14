function plotPSAsplit(splitValues, metricVector, trl, mat, figTitle, plotTitle, printCount)
    if nargin < 7, printCount = true; end
    edges = sort(splitValues(:)'); %#ok<TRSRT> 
    numSplits = length(edges);
    numBins = numSplits + 1;
    
    % Aloca matrizes separadas para Acurácia, Sensibilidade e Contagens
    barMatrixAcc  = zeros(numBins, 2);
    barMatrixSens = zeros(numBins, 2);
    countMatrix   = zeros(numBins, 2);
    xlabelPSA     = cell(1, numBins);
    
    calcDPrime = @(hit, total) norminv(min(max(hit./total, 0.001), 0.999)) * sqrt(2);
    
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
        
        % Se houver dados, calcula counts, porcentagens e d-primes
        if any(mask)
            [~, countsSub] = getPSAeffect(trl(mask));
            
            % Salva os denominadores reais [Sacádico, Não-sacádico]
            countMatrix(bIdx, :) = [countsSub(2,2), countsSub(2,3)];
            
            % Acurácia (%)
            pct_s = (countsSub(1,2) / countsSub(2,2)) * 100;
            pct_n = (countsSub(1,3) / countsSub(2,3)) * 100;
            barMatrixAcc(bIdx, :) = [pct_s, pct_n];
            
            % Sensibilidade (d')
            d_s = calcDPrime(countsSub(1,2), countsSub(2,2));
            d_n = calcDPrime(countsSub(1,3), countsSub(2,3));
            barMatrixSens(bIdx, :) = [d_s, d_n];
        else
            barMatrixAcc(bIdx, :)  = [0, 0];
            barMatrixSens(bIdx, :) = [0, 0];
            countMatrix(bIdx, :)   = [0, 0];
        end
    end
    
    if ~printCount, countMatrix = []; end
    renderPSAsplit(barMatrixAcc, countMatrix, numBins, xlabelPSA, 'Acertos (%)', ...
        [figTitle, ' - Accuracy'], [plotTitle, ' (Acurácia)'], mat, false);
    renderPSAsplit(barMatrixSens, countMatrix, numBins, xlabelPSA, 'Sensibilidade (d'')', ...
        [figTitle, ' - Sensitivity'], [plotTitle, ' (Sensibilidade)'], mat, true);
end
function plotPSAsplit(splitValues, metricVector, trl, mat, figTitle, plotTitle, printCount)
    if nargin < 7, printCount = true; end
    
    % --- ZERO-AWARE QUANTILE HANDLING ---
    hasZeroInflation = (mean(metricVector == 0) >= 0.20);
    
    if hasZeroInflation
        posMetrics = metricVector(metricVector > 0);
        numSplits  = length(splitValues); 
        
        % Calculate quantiles ONLY on positive data
        pSteps = linspace(0, 1, numSplits + 2); 
        posEdges = quantile(posMetrics, pSteps(2:end-1));
        posEdges = unique(posEdges);
        
        % 1 bin for zeros + bins for positive data
        numBins = length(posEdges) + 2; 
    else
        edges   = sort(splitValues(:)');
        edges   = unique(edges);
        numBins = length(edges) + 1;
    end
    
    % Preallocate matrices
    barMatrixAcc  = zeros(numBins, 2);
    barMatrixSens = zeros(numBins, 2);
    countMatrix   = zeros(numBins, 2);
    xlabelPSA     = cell(1, numBins);
    
    calcDPrime = @(hit, total) norminv(min(max(hit./total, 0.001), 0.999)) * sqrt(2);
    
    for bIdx = 1:numBins
        if hasZeroInflation
            % --- ZERO-INFLATED MASK LOGIC ---
            if bIdx == 1
                % Bin 1: Strictly Zero
                mask = (metricVector == 0);
                xlabelPSA{bIdx} = '= 0';
                
            elseif bIdx == 2
                % Bin 2: First positive bin (strictly > 0 to avoid double-counting)
                mask = (metricVector > 0 & metricVector < posEdges(1));
                xlabelPSA{bIdx} = sprintf('(0, %.3f)', posEdges(1));
                
            elseif bIdx == numBins
                % Last Bin: Top positive quantile
                mask = (metricVector >= posEdges(end));
                xlabelPSA{bIdx} = sprintf('\\geq %.3f', posEdges(end));
                
            else
                % Intermediate positive bins
                eIdx = bIdx - 2; % Shift index for posEdges
                mask = (metricVector >= posEdges(eIdx) & metricVector < posEdges(eIdx+1));
                xlabelPSA{bIdx} = sprintf('[%.3f, %.3f)', posEdges(eIdx), posEdges(eIdx+1));
            end
            
        else
            % --- STANDARD CONTINUOUS MASK LOGIC ---
            if bIdx == 1
                mask = (metricVector < edges(1));
                xlabelPSA{bIdx} = sprintf('< %.3f', edges(1));
                
            elseif bIdx == numBins
                mask = (metricVector >= edges(end));
                xlabelPSA{bIdx} = sprintf('\\geq %.3f', edges(end));
                
            else
                mask = (metricVector >= edges(bIdx-1) & metricVector < edges(bIdx));
                xlabelPSA{bIdx} = sprintf('[%.3f, %.3f)', edges(bIdx-1), edges(bIdx));
            end
        end
        
        % --- Calculate Statistics ---
        if any(mask)
            [~, countsSub] = getPSAeffect(trl(mask));
            countMatrix(bIdx, :) = [countsSub(2,2), countsSub(2,3)];
            
            pct_s = (countsSub(1,2) / countsSub(2,2)) * 100;
            pct_n = (countsSub(1,3) / countsSub(2,3)) * 100;
            barMatrixAcc(bIdx, :) = [pct_s, pct_n];
            
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
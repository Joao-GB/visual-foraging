function [dPrimeMap, countMap, X, Y] = computeSpatialDPrime(trl, coordsSacc, coordsNSacc, gridLim, gridRes, radius, minTrials)
    % Cria a malha de coordenadas
    xVec = -gridLim:gridRes:gridLim;
    yVec = -gridLim:gridRes:gridLim;
    [X, Y] = meshgrid(xVec, yVec);
    
    dPrimeMap = nan(size(X));
    countMap  = zeros(size(X));
    
    calcDPrime = @(hit, total) norminv(min(max(hit/total, 0.001), 0.999)) * sqrt(2);
    
    % Loop sobre cada "pixel" do grid
    for r = 1:size(X, 1)
        for c = 1:size(X, 2)
            px = X(r, c);
            py = Y(r, c);
            
            % Distância euclidiana de todos os probes para este pixel
            distSacc  = sqrt((coordsSacc(:,1) - px).^2 + (coordsSacc(:,2) - py).^2);
            distNSacc = sqrt((coordsNSacc(:,1) - px).^2 + (coordsNSacc(:,2) - py).^2);
            
            % Máscaras lógicas de quais trials caem dentro do raio
            maskSacc  = distSacc <= radius;
            maskNSacc = distNSacc <= radius;
            
            hitsTotal = 0;
            trialsTotal = 0;
            
            % Truque: usamos sua própria função getPSAeffect como "caixa preta"
            % 1. Extrair acertos Sacádicos locais
            if any(maskSacc)
                [~, cSacc] = getPSAeffect(trl(maskSacc));
                hitsTotal   = hitsTotal + cSacc(1,2);
                trialsTotal = trialsTotal + cSacc(2,2);
            end
            
            % 2. Extrair acertos Não-sacádicos locais
            if any(maskNSacc)
                [~, cNSacc] = getPSAeffect(trl(maskNSacc));
                hitsTotal   = hitsTotal + cNSacc(1,3);
                trialsTotal = trialsTotal + cNSacc(2,3);
            end
            
            % 3. Calcular o d' local se tivermos trials suficientes
            countMap(r, c) = trialsTotal;
            if trialsTotal >= minTrials
                dPrimeMap(r, c) = calcDPrime(hitsTotal, trialsTotal);
            end
        end
    end
    
    % Encontra os índices dos pixels que passaram do limite mínimo
    validIdx = ~isnan(dPrimeMap);
    
    if any(validIdx, 'all')
        % Cria os interpoladores ('none' impede que ele extrapole para fora da área convexa)
        F_dprime = scatteredInterpolant(X(validIdx), Y(validIdx), dPrimeMap(validIdx), 'linear', 'none');
        F_counts = scatteredInterpolant(X(validIdx), Y(validIdx), countMap(validIdx), 'linear', 'none');
        
        % Preenche os buracos internos mantendo a forma convexa
        dPrimeMap = F_dprime(X, Y);
        countMap  = F_counts(X, Y);
    end
end
function plotPSAdprimeMap(pPos, pFix, nPos, trl, ~, rotated)
    if nargin < 6, rotated = false; end
    
    % Parâmetros do Mapa de Calor
    gridLim = 10;      % Limite do mapa em graus visuais (-12 a 12)
    gridRes = .5;     % Resolução do grid (tamanho do pixel em DVA)
    radius  = 2.5;     % Raio da janela deslizante para agrupamento de trials
    minTrials = 10;    % Mínimo de trials numa região para calcular o d' (evita ruído)
    
    % 2. MAPA 1: Relativo à Posição Real do Alvo (probePos = 0,0)
    % Para este mapa, o alvo sacádico está SEMPRE em (0,0).
    relSacc1 = zeros(size(pPos)); % (pPos - pPos)
    relNSacc1 = nPos - pPos;
    
    [mapDPrime1, mapCounts1, X, Y] = computeSpatialDPrime(trl, relSacc1, relNSacc1, ...
        gridLim, gridRes, radius, minTrials);
        
    % 3. MAPA 2: Relativo à Posição Final do Olho (probePosFix = 0,0)
    relSacc2 = pPos - pFix;
    relNSacc2 = nPos - pFix;
    
    [mapDPrime2, mapCounts2, ~, ~] = computeSpatialDPrime(trl, relSacc2, relNSacc2, ...
        gridLim, gridRes, radius, minTrials);

    % 4. Renderização
    if rotated, mapTitle = 'PSA d-prime saccade-aligned spatial map';
    else,       mapTitle = 'PSA d-prime spatial map';
    end
   renderPSAdprimeMaps(X, Y, mapDPrime1, mapCounts1, mapDPrime2, mapCounts2, ...
        '(0,0) = Centro do probe', ...
        '(0,0) = Fixação no probe', ...
        mapTitle);

end
function renderPSAdprimeMaps(X, Y, dPrime1, count1, dPrime2, count2, sub1, sub2, figName)
    f = figure('Name', figName, 'Color', 'w', 'Position', [100 100 1100 550]);
    
    % Encontra os limites globais do d-prime para unificar a escala de cores
    minD = min([dPrime1(:); dPrime2(:)], [], 'omitnan');
    maxD = max([dPrime1(:); dPrime2(:)], [], 'omitnan');
    if isnan(minD) || isnan(maxD), minD = 0; maxD = 1; end
    if minD == maxD, maxD = minD + 1; end % Evita erro no eixo de cores
    
    % Cria o layout lado a lado
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    nexttile;
    drawSingleMap(X, Y, dPrime1, count1, sub1, [minD, maxD]);
    
    % --- Subplot 2 ---
    ax2 = nexttile;
    drawSingleMap(X, Y, dPrime2, count2, sub2, [minD, maxD]);
    
    % Adiciona UMA colorbar compartilhada na direita do layout
    cb = colorbar;
    cb.Layout.Tile = 'east'; 
    ylabel(cb, 'Sensibilidade (d'')', 'FontSize', 12, 'FontWeight', 'bold');
    
    colormap(f, parula); % Aplica o mapa de cores à figura inteira
    
    % Título Geral
%     title(t, 'Distribuição Espacial de Sensibilidade (d'')', 'FontSize', 16, 'FontWeight', 'bold');
end

% --- Função auxiliar para não repetir o código de renderização ---
function drawSingleMap(X, Y, dPrimeMap, countMap, subTitle, cLimits)
    % Renderiza o mapa de calor
    h = imagesc(X(1,:), Y(:,1)', dPrimeMap);
    
    % Escala de transparência (Alpha)
    maxReliableTrials = 30;
    alphaMap = countMap / maxReliableTrials;
    alphaMap(alphaMap > 1) = 1; 
    alphaMap(isnan(dPrimeMap)) = 0; 
    set(h, 'AlphaData', alphaMap);
    
    set(gca, 'YDir', 'normal'); 
    
    % Unifica os limites de cor (clim é usado no MATLAB R2022a+, ou use caxis em versões mais antigas)
    clim(cLimits); 
    
    hold on;
    [~, hc] = contour(X, Y, dPrimeMap, 5, 'LineColor', [1 1 1], 'LineWidth', 0.5, 'LineStyle', '-');
    plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    
    % Estética
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', [0.95 0.95 0.95]);
    grid on;
    axis equal tight;
    
    xlabel('Distância Horizontal (dva)');
    ylabel('Distância Vertical (dva)');
    title(subTitle, 'FontSize', 13);
end
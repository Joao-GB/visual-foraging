function plotPSASaccProps(preProbePos, probePos, preProbePosFix, probePosFix, drP)
    % plotPSASaccProps para métricas de execução de sacadas

    vActualSacc = probePosFix - preProbePosFix;
%     vTargetSacc = probePos - preProbePos;
    
    % Direção e amplitude da sacada
    actualSaccAngle = angle(complex(vActualSacc(:,1), vActualSacc(:,2)));
    actualSaccAmp   = vecnorm(vActualSacc, 2, 2);
%     targetSaccAmp   = vecnorm(vTargetSacc, 2, 2); % Ideal amplitude for reference
    
    % Erro da sacada
    saccError = vecnorm(probePos - probePosFix, 2, 2);

    figure('Name', 'PSA Sacc Props', 'Color', 'w', 'Position', [50 50 700 700]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    alphaLevel = 0.5;
    colorPostFix = drP.darkGreen;
    colorAmp     = drP.darkGreen;
    colorError   = repmat(drP.whiteGrey, [1 3]); 
    
    % Subplot 1: Histograma polar da direção sacádica
    nexttile;
    numBins = 20;
    binEdges = getCardinalCenteredEdges(numBins);
    
    polarhistogram(actualSaccAngle, binEdges, 'FaceColor', colorPostFix, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    title('Direção das sacadas');
    set(gca, 'ThetaZeroLocation', 'right', 'ThetaDir', 'counterclockwise'); 

    % Subplot 2: Amplitude e direção sacádicas
    nexttile;
    polarscatter(actualSaccAngle, actualSaccAmp, 25, colorPostFix, 'filled', 'MarkerFaceAlpha', alphaLevel);
    title('Direção e amplitude das sacadas');
    set(gca, 'ThetaZeroLocation', 'right', 'ThetaDir', 'counterclockwise'); 
    rLim = rlim;
    rlim([0 rLim(2)]);

    % Subplot 3: amplitude das sacadas
    nexttile;
    hold on;
    histogram(actualSaccAmp, 20, 'FaceColor', colorAmp, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    
%     meanTarget = mean(targetSaccAmp);
%     xline(meanTarget, 'r--', 'LineWidth', 1.5, 'Label', 'Target Mean');
    
    xlabel('Amplitude (dva)');
    ylabel('Contagem de trials');
    title({'Distribuição da amplitude das sacadas'});
    grid on;
    set(gca, 'TickDir', 'out', 'Box', 'off');
    yl = ylim;
    ylim([-0.003*diff(yl), yl(2)])

    % Subplot 4: distribuição do erro de sacada
    nexttile;
    histogram(saccError, 20, 'FaceColor', colorError, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    xlabel('Erro (dva)');
    ylabel('Contagem de trials');
    title({'Erro da sacada em relação ao probe'});
    grid on;
    set(gca, 'TickDir', 'out', 'Box', 'off');

    yl = ylim;
    ylim([-0.003*diff(yl), yl(2)])
end

function edges = getCardinalCenteredEdges(numBins)
    % O número de bins deve ser múltiplo de 4
    
    binWidthDeg = 360 / numBins;
    halfBinDeg  = binWidthDeg / 2;
    
    edgesDeg = -halfBinDeg : binWidthDeg : (360 - halfBinDeg);
    edges = deg2rad(edgesDeg);
end
function foragingEyePlots(fixPos, fixPosFix, probePos, probePosFix, nSaccProbePos)
    % foragingEyePlot para vizualizar a geometria do comportamento oculomotor
    % Inputs devem ser matrizes Nx2 em dva
    
    % 1. Todas posições são realinhadas de acordo com o novo (0,0), o centro
    %    do estímulo pré-probe
    vProbePos    = probePos - fixPos;        % Posição do probe sacádico
    vNSaccPos    = nSaccProbePos - fixPos;   % Posição do probe não-sacádico
    vPreProbeFix = fixPosFix - fixPos;       % Fixação do olho em pré-probe
    vProbeFix    = probePosFix - fixPos;     % Fixação do olho em probe
    
    % 2. Representação complexa dos vetores relevantes
    zProbePos = complex(vProbePos(:,1), vProbePos(:,2));
    zNSaccPos = complex(vNSaccPos(:,1), vNSaccPos(:,2));
    zProbeFix = complex(vProbeFix(:,1), vProbeFix(:,2));
    
    % Rotação sem alterar a amplitude (distância em dva)
    % Dividimos apenas pela "direção" (vetor unitário) do probe sacádico
    uProbePos  = exp(1i * angle(zProbePos)); 
    rotNSacc   = zNSaccPos ./ uProbePos;
    rotLanding = zProbeFix ./ uProbePos;
    rotProbe   = zProbePos ./ uProbePos; % Fica sempre no eixo x, mas mantém a distância original
    
    % 3. Ângulo do movimento ocular (entre fixação probe e fixação pré-probe)
    vActualSacc = vProbeFix - vPreProbeFix;
    actualSaccAngle = angle(complex(vActualSacc(:,1), vActualSacc(:,2)));
    
    % 4. Figura em si
    figure('Name', 'foragingEyePlots', 'Position', [100 200 1200 400]);
    t = tiledlayout(1,3);
    alphaLevel = 0.5;
    
    colorPreFix  = [0.3 0.3 0.3];    % cinza escuro
    colorPostFix = [0.2 0.8 0.2];    % verde 
    colorProbe   = [0.0 0.45 0.74];  % azul
    colorNSacc   = [0.85 0.33 0.10]; % laranja
    
    %% Subplot 1: espaço inteiro em dva
    nexttile; hold on;
    
    % Linhas conectando probePos a probeFix (para ter noção do erro da
    % sacada à probe)
    x_lines = [vProbePos(:,1)'; vProbeFix(:,1)'; NaN(1, size(vProbePos,1))];
    y_lines = [vProbePos(:,2)'; vProbeFix(:,2)'; NaN(1, size(vProbePos,1))];
    plot(x_lines(:), y_lines(:), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, 'HandleVisibility', 'off');

    % Posição dos estímulos
    h2 = scatter(vProbePos(:,1), vProbePos(:,2), 35, colorProbe, 'filled', 'MarkerFaceAlpha', alphaLevel, 'MarkerEdgeColor', 'flat', 'LineWidth', 1);
    h3 = scatter(vNSaccPos(:,1), vNSaccPos(:,2), 35, colorNSacc, 'filled', 'MarkerFaceAlpha', alphaLevel, 'MarkerEdgeColor', 'flat', 'LineWidth', 1);
    
    % Posição ocular
    h4 = scatter(vPreProbeFix(:,1), vPreProbeFix(:,2), 25, colorPreFix, 'filled', 'MarkerFaceAlpha', alphaLevel);
    h5 = scatter(vProbeFix(:,1), vProbeFix(:,2), 25, colorPostFix, 'filled', 'MarkerFaceAlpha', alphaLevel);

    % Origem
    h1 = plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    
    axis equal; grid on; set(gca, 'TickDir', 'out', 'Box', 'off');
    xlabel('Horizontal (dva)'); ylabel('Vertical (dva)');
    title({'Posição de estímulos e fixações', 'em relação à pré-probe'});
    
    %% Subplot 2: diração das sacadas em direção ao probe
    nexttile;
    polarhistogram(actualSaccAngle, 18, 'FaceColor', colorPostFix, 'FaceAlpha', 0.7);
    title({'Direção das sacadas aos probes'});
    set(gca, 'ThetaZeroLocation', 'right', 'ThetaDir', 'counterclockwise');
    
    %% Subplot 3: posição relativa do probe não sacádico
    nexttile; hold on;

    % Linhas entre rotNSacc e rotLanding no espaço rotacionado
    x_lines3 = [real(rotNSacc)'; real(rotLanding)'; NaN(1, length(rotNSacc))];
    y_lines3 = [imag(rotNSacc)'; imag(rotLanding)'; NaN(1, length(rotNSacc))];
    plot(x_lines3(:), y_lines3(:), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, 'HandleVisibility', 'off');

    plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    scatter(real(rotProbe), imag(rotProbe), 35, colorProbe, 'filled', 'MarkerFaceAlpha', alphaLevel, 'MarkerEdgeColor', 'flat', 'LineWidth', 1);
    scatter(real(rotNSacc), imag(rotNSacc), 35, colorNSacc, 'filled', 'MarkerFaceAlpha', alphaLevel, 'MarkerEdgeColor', 'flat', 'LineWidth', 1);
    scatter(real(rotLanding), imag(rotLanding), 25, colorPostFix, 'filled', 'MarkerFaceAlpha', alphaLevel);
    
    
    axis equal; grid on; set(gca, 'TickDir', 'out', 'Box', 'off');
    xlabel('Horizontal transformado (dva)'); 
    ylabel('Horizontal transformado (dva)');
    title({'Posição de probes não-sacádicos', 'em relação à direção das sacadas'});

    lgd = legend([h1 h2 h3 h4 h5], {'pré-probe-pos', 'probe pos', 'nSacc probe pos', 'pré-probe fix', 'probe fix'});
    lgd.Layout.Tile = 'south'; 
    lgd.Orientation = 'horizontal';
end
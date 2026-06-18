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
    rotProbe   = zProbePos ./ uProbePos; % Estará sempre no eixo X (mas mantendo a distância original)
    
    % --- 3. POLAR ANGLES ---
    % Ângulo do movimento ocular real (Fixação Probe em relação à Fixação Pré-Probe)
    vActualSacc = vProbeFix - vPreProbeFix;
    actualSaccAngle = angle(complex(vActualSacc(:,1), vActualSacc(:,2)));
    
    % --- 4. VISUALIZATION ---
    figure('Name', 'Geometria Espacial e Oculomotora', 'NumberTitle', 'off', 'Position', [100 200 1200 400]);
    alphaLevel = 0.5; % Aumentado para melhor visibilidade dos pontos
    
    % Definindo cores mais distintas
    colorPreFix  = [0.2 0.2 0.2];    % Cinza escuro (quase preto)
    colorPostFix = [0.2 0.8 0.2];    % Verde vibrante
    colorProbe   = [0.0 0.45 0.74];  % Azul
    colorNSacc   = [0.85 0.33 0.10]; % Laranja/Vermelho
    
    % --- Subplot 1: Absolute Retinal Space ---
    subplot(1,3,1); hold on;
    
    % Linhas conectando probe pos a probe fix (Trajetória da Sacada)
    % Truque de vetorização com NaN para desenhar todas as linhas sem usar loop
    x_lines = [vProbePos(:,1)'; vProbeFix(:,1)'; NaN(1, size(vProbePos,1))];
    y_lines = [vProbePos(:,2)'; vProbeFix(:,2)'; NaN(1, size(vProbePos,1))];
    plot(x_lines(:), y_lines(:), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    % Eye positions (Behavior)
    scatter(vPreProbeFix(:,1), vPreProbeFix(:,2), 25, colorPreFix, 'filled', 'MarkerFaceAlpha', alphaLevel);
    scatter(vProbeFix(:,1), vProbeFix(:,2), 25, colorPostFix, 'filled', 'MarkerFaceAlpha', alphaLevel);
    % Stimulus positions (Geometry)
    scatter(vProbePos(:,1), vProbePos(:,2), 35, colorProbe, 'filled', 'MarkerFaceAlpha', alphaLevel);
    scatter(vNSaccPos(:,1), vNSaccPos(:,2), 35, colorNSacc, 'filled', 'MarkerFaceAlpha', alphaLevel);
    
    % Origin marker
    plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    
    axis equal; grid on; set(gca, 'TickDir', 'out', 'Box', 'off');
    xlabel('Horizontal (dva)'); ylabel('Vertical (dva)');
    title('Coordenadas Absolutas Centradas');
    legend('Pre-probe fix', 'Probe fix', 'Probe pos', 'N-sacc probe pos', 'Origem (0,0)', 'Location', 'best');
    
    % --- Subplot 2: Polar Direction Distribution ---
    subplot(1,3,2);
    % Histograma polar da sacada real
    polarhistogram(actualSaccAngle, 18, 'FaceColor', colorPostFix, 'FaceAlpha', 0.7);
    title({'Direção da Sacada', '(Probe Fix rel. Pre-Probe Fix)'});
    set(gca, 'ThetaZeroLocation', 'right', 'ThetaDir', 'counterclockwise');
    
    % --- Subplot 3: Saccade-Rotated Space ---
    subplot(1,3,3); hold on;
    % Rotated eye landings
    scatter(real(rotLanding), imag(rotLanding), 25, colorPostFix, 'filled', 'MarkerFaceAlpha', alphaLevel);
    % Rotated Non-Saccadic probe
    scatter(real(rotNSacc), imag(rotNSacc), 35, colorNSacc, 'filled', 'MarkerFaceAlpha', alphaLevel);
    % Rotated Saccadic probe (Variação de amplitude real no eixo X)
    scatter(real(rotProbe), imag(rotProbe), 35, colorProbe, 'filled', 'MarkerFaceAlpha', alphaLevel);
    
    % Reference markers: Origin (Fixation)
    plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    
    axis equal; grid on; set(gca, 'TickDir', 'out', 'Box', 'off');
    xlabel('X Rotacionado (dva)'); 
    ylabel('Y Rotacionado (dva)');
    title('Espaço Rotacionado pela Saccade');
end
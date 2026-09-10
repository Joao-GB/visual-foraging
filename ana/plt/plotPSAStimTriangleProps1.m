function [group1, group2] = plotPSAStimTriangleProps1(preProbePos, probePos, nSaccProbePos, drP, hypothesis, doPlot)

% plotPSAStimTriangleProps1
%
% Me interessam apenas posições, e não fixações, para essa função.
%
% As extremidades do triângulo são denotadas por:
%
%   P = pré-probe / fixation
%   S = probe sacádico
%   N = probe não-sacádico
%
% O espaço geométrico é definido por três variáveis:
%
%   PN    = distância entre o pré-probe e o probe não-sacádico
%   PS    = distância entre o pré-probe e o probe sacádico
%   theta = ângulo NPS, formado no pré-probe
%
% hypothesis:
%
%   0 = no hypothesis-specific coloring
%   1 = direction hypothesis
%   2 = eccentricity hypothesis
%   3 = spotlight hypothesis
%
% A figura mostra as três projeções 2D desse espaço 3D:
%
%   1) PN x theta
%   2) PS x theta
%   3) PN x PS
%
% E, na posição (2,2), uma visualização 3D:
%
%   x = PN
%   y = PS
%   z = theta

%% Default hypothesis

if nargin < 5, hypothesis = 0; end
if nargin < 6, doPlot = true;  end

%% Cálculo das distâncias euclidianas

PS = vecnorm(probePos - preProbePos, 2, 2);
PN = vecnorm(nSaccProbePos - preProbePos, 2, 2);

% Distância entre os dois probes
NS = vecnorm(probePos - nSaccProbePos, 2, 2);

% Garantir vetores coluna
PS = PS(:);
PN = PN(:);
NS = NS(:);

%% Ângulo formado no pré-probe

theta = acosd( ...
    (PN.^2 + PS.^2 - NS.^2) ./ (2 .* PN .* PS) ...
);

%% Métricas auxiliares

% Relação entre as duas excentricidades.
%
%   1   = exatamente iguais
%   1.2 = diferença de até 20%
%   1.5 = diferença de 50%
%
eccRatio = max(PN, PS) ./ min(PN, PS);

%% Definição dos grupos de acordo com a hipótese

% group1 / group2 são os dois grupos que serão destacados.
% selected indica quais trials pertencem a algum dos grupos.

group1 = false(size(PN));
group2 = false(size(PN));

switch hypothesis

    case 0
        % Nenhum destaque específico

    case 1
        % -------------------------------------------------------------
        % H1: DIRECTION
        %
        % Primeiro restringimos aos trials em que PN e PS são
        % aproximadamente iguais. Assim, theta realmente representa
        % a diferença de direção entre os probes, sem confundir isso
        % com diferenças de excentricidade.
        %
        % Similar eccentricity:
        %       PS/PN dentro de aproximadamente +/- 20%
        %
        % Same direction:
        %       theta < 20 deg
        %
        % Opposite direction:
        %       theta > 160 deg
        % -------------------------------------------------------------

        similarEcc = eccRatio <= 2;

        group1 = similarEcc & theta < 40;
        group2 = similarEcc & theta > 140;

    case 2
        % -------------------------------------------------------------
        % H2: ECCENTRICITY
        %
        % Comparamos trials com eccentricidades semelhantes e
        % claramente diferentes.
        %
        % Similar:
        %       eccRatio <= 1.2
        %
        % Dissimilar:
        %       eccRatio >= 1.5
        %
        % A região intermediária é deixada sem destaque.
        % -------------------------------------------------------------

        group1 = eccRatio <= 1.25;
        group2 = eccRatio >= 1.75;

    case 3
        % -------------------------------------------------------------
        % H3: SPOTLIGHT
        %
        % Comparamos trials em que os probes estão próximos ou
        % distantes entre si.
        %
        % Close:
        %       NS < 5 dva
        %
        % Far:
        %       NS >= 5 dva
        % -------------------------------------------------------------

        H3cut = 8;
        group1 = NS < H3cut;
        group2 = NS >= H3cut;

    otherwise
        error('hypothesis must be 0, 1, 2, or 3.');

end

%% Cores

alphaLevel = 0.5;
colorMetric = drP.orange;

% Cores usadas apenas para o destaque das hipóteses
colorGroup1 = [0.8500 0.3250 0.0980];
colorGroup2 = [0.0000 0.4470 0.7410];
colorOther  = [0.75 0.75 0.75];

%% Figura
if doPlot
    figure( ...
        'Name', 'PSA Triangle Props 1', ...
        'Color', 'w', ...
        'Position', [100 100 1300 750] ...
    );
    
    tiledlayout(2, 3, 'TileSpacing', 'compact');
    
    %% Subplot 1: PN x theta
    
    nexttile
    
    hold on;
    
    scatter( ...
        theta, PN, ...
        25, colorMetric, ...
        'filled', ...
        'MarkerFaceAlpha', alphaLevel ...
    );
    
    xlabel('Ângulo formado no pré-probe (deg)');
    ylabel('Distância P-N (dva)');
    title('Ângulo × probe não-sacádico');
    
    grid on;
    set(gca, 'TickDir', 'out', 'Box', 'off');
    
    xlim([0 180]);
    
    %% Subplot 2: PS x theta
    
    nexttile
    
    hold on;
    
    scatter( ...
        theta, PS, ...
        25, colorMetric, ...
        'filled', ...
        'MarkerFaceAlpha', alphaLevel ...
    );
    
    xlabel('Ângulo formado no pré-probe (deg)');
    ylabel('Distância P-S (dva)');
    title('Ângulo × probe sacádico');
    
    grid on;
    set(gca, 'TickDir', 'out', 'Box', 'off');
    
    xlim([0 180]);
    
    %% Subplot 3: PN x PS
    
    nexttile
    
    hold on;
    
    scatter( ...
        PN, PS, ...
        25, colorMetric, ...
        'filled', ...
        'MarkerFaceAlpha', alphaLevel ...
    );
    
    % Determinar limites comuns para manter a diagonal y = x
    
    maxVal = max([PN; PS]);
    
    xlim([0 maxVal]);
    ylim([0 maxVal]);
    
    % Linha de igualdade PN = PS
    
    plot( ...
        [0 maxVal], [0 maxVal], ...
        'k--', ...
        'LineWidth', 1 ...
    );
    
    xlabel('Distância P-N (dva)');
    ylabel('Distância P-S (dva)');
    title('Probe não-sacádico × probe sacádico');
    
    axis equal;
    grid on;
    set(gca, 'TickDir', 'out', 'Box', 'off');
    
    %% Subplot 4: 3D space
    
    nexttile(5)
    
    hold on;
    
    % -------------------------------------------------------------
    % Primeiro, se houver uma hipótese selecionada, mostrar todos
    % os trials em cinza claro.
    % -------------------------------------------------------------
    
    if hypothesis == 0
    
        % Sem hipótese: comportamento original
        scatter3( ...
            PN, PS, theta, ...
            25, colorMetric, ...
            'filled', ...
            'MarkerFaceAlpha', alphaLevel ...
        );
    
    else
    
        % Trials não selecionados
        other = ~(group1 | group2);
    
        scatter3( ...
            PN(other), PS(other), theta(other), ...
            20, colorOther, ...
            'filled', ...
            'MarkerFaceAlpha', 0.25 ...
        );
    
        % Grupo 1
        scatter3( ...
            PN(group1), PS(group1), theta(group1), ...
            35, colorGroup1, ...
            'filled', ...
            'MarkerFaceAlpha', 0.8 ...
        );
    
        % Grupo 2
        scatter3( ...
            PN(group2), PS(group2), theta(group2), ...
            35, colorGroup2, ...
            'filled', ...
            'MarkerFaceAlpha', 0.8 ...
        );
    
    end
    
    xlabel('Distância P-N (dva)');
    ylabel('Distância P-S (dva)');
    zlabel('Ângulo no pré-probe (deg)');
    
    switch hypothesis
    
        case 0
            title('Espaço 3D');
    
        case 1
            title({'Espaço 3D', 'Direção: mesmo vs. oposto'});
    
        case 2
            title({'Espaço 3D', 'Excentricidade: similar vs. diferente'});
    
        case 3
            title("Espaço 3D - Spotlight: NS < " + H3cut + ...
          " vs. NS \geq " + H3cut + " dva");
    
    end
    
    grid on;
    set(gca, 'TickDir', 'out', 'Box', 'off');
    
    view(3);
end
end
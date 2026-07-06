function plotPSAStimTriangleProps(preProbePos, probePos, nSaccProbePos, drP)
% Me interessam apenas posições, e não fixações, para essa função
% The function produces two standard scatter plots using the color code from drP:
%   1) Shape coordinates:
%        x = angle opposite NS
%        y = (PN - PS)/NS
%
%   2) Scale-invariant side ratios:
%        x = PN/NS
%        y = PS/NS

% Cálculo das distâncias euclidianas (lados do triângulo)
PS = vecnorm(probePos - preProbePos, 2, 2);
PN = vecnorm(nSaccProbePos - preProbePos, 2, 2);
NS = vecnorm(probePos - nSaccProbePos, 2, 2);

NS = NS(:);
PN = PN(:);
PS = PS(:);

%% Métricas apresentadas na figura
shape1 = acosd( (PN.^2 + PS.^2 - NS.^2) ./ (2*PN.*PS) );
shape2 = (PS - PN)./NS;
ratio1 = PN./NS;
ratio2 = PS./NS;

figure('Name', 'PSA Triangle Props', 'Color','w', 'Position', [100 200 1300 400]);
tiledlayout(1,3,'TileSpacing','compact');

alphaLevel = 0.5;
colorMetric = drP.orange;

%% Subplot 1: ângulo por diferença entre os lados
% Entende-se como o quão aberto e o quão assimétrico é o triângulo formado
% por eles
% Interpretação:
% - Muitos pontos ao longo de zero: triângulos isósceles com vários ângulos;
% - Ângulo nulo: colineares;
% - Pontos concentrados em -1: probe sacádico na metade do caminho em
%                               relação ao não-sacádico;
% - Pontos concentrados em 1:  probe não sacádico 2 vezes mais perto que o
%                               sacádico

nexttile
hold on;
scatter(shape1, shape2, 25, colorMetric, 'filled', 'MarkerFaceAlpha', alphaLevel);

xlabel('Ângulo formado no pré-probe (deg)')
ylabel('Assimetria dos vetores de probes')
title('Coordenadas de formato')
grid on
set(gca, 'TickDir', 'out', 'Box', 'off');
xlim([0 180]);

%% Subplot 2
% Interpretação:
% - Na diagonal: ambos os vetores com mesmo tamanho;
% - Posição ao longo da diagonal: perto de 0 indica probes muito distantes
%                                  entre si; longe indica que ambos os
%                                  vetores são bem maiores que NS. Logo,
%                                  abaixo de raio 1 para ângulos obtusos,
%                                  acima para ângulos agudos.
% - Abaixo da diagonal: X > Y logo vetor não-sacádico maior que o sacádico;
% - Acima da diagonal:  X < Y logo vetor não-sacádico menor que o sacádico;
% -

nexttile
hold on;
scatter(ratio1, ratio2, 25, colorMetric, 'filled', 'MarkerFaceAlpha', alphaLevel);

xlabel('Vetor não-sacádico')
ylabel('Vetor sacádico')
title({'Comprimento dos vetores normalizados', 'pela distância entre suas pontas'})
axis equal
grid on
set(gca, 'TickDir', 'out', 'Box', 'off');

maxVal = max([xlim, ylim]); 
xlim([0 maxVal]);
ylim([0 maxVal]);

nexttile
histogram(NS, 20, 'FaceColor', colorMetric, 'EdgeColor', 'w', 'FaceAlpha', alphaLevel);
xlabel('Distância (dva)')
ylabel('Contagem de trials')
title('Distância absoluta entre probes')
grid on
set(gca, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top');

end
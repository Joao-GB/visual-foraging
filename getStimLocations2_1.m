function [fixCenter, stimCenters, rMax] = getStimLocations2_1(screenRes, ROIparams, nStims, minDist, rStim, c)
% Poisson disk sampling com raio de busca de vizinhos restrito ao 
% intervalo [minDist, minDist*(1+c)]. 
% A nova lógica: 
% (1) Preenche o retÂngulo dado por [0, 0, screenRes(1), screenRes(2)] 
%     com o ruído azul que respeita os params minDist e c; 
% (2) ROIparams e 1x3 (ou 1x4) com formato (x, y, 0 ou 1, r), em que x e 1 
%     são as coordenadas do centro do ROI e 0 para ROI circular, 1 para 
%     quadrado; r, opcional, indica o raio do círculo (ou metade do lado 
%     do quadrado) 
% (3) nStims e um inteiro que indica a quantidade de estimulos que 
%     serão incluídos no ROI. Se o raio não for dado, o ROI devera ser 
%     expandido até que inclua nStims pontos; se for dado, devem ser 
%     incluídos apenas os que estejam dentro do ROI e distem ao menos 
%     minDist de suas bordas; se nStims = [] e o raio tambem for vazio, 
%     o raio do ROI passa a ser o maior quadrado (ou círculo) que cabe 
%     dentro da tela. 
% (4) Ao todo, devem ser exatamente nStims pontos em stimCenters e 
%     fixCenter de ser um dos stimCenters escolhidos ao acaso
%
% Para testar:
% ROIparams = [960 540 1]; nStims=10;
% mW = 52.4; sR = 1920; sD = 57;
% minDist_dva = 4.5; rStim_dva = 1.5/2;
% minDist = dva2pix(sD, mW, sR, minDist_dva); rStim = dva2pix(sD, mW, sR, rStim_dva);
% screenRes = [1920 1080];
% [fixCenter, stimCenters] = getStimLocations2_1(screenRes, ROIparams, nStims, minDist, rStim, 1);
% figure;scatter(stimCenters(1,:), stimCenters(2,:)); hold on;rectangle('Position', [0 0 screenRes(1) screenRes(2)], 'EdgeColor', 'r', 'LineWidth', 2);axis equal;

if nargin < 5, rStim = minDist/2; end
if nargin < 6, c = 1; end

% Extrai os parâmetros
cx = ROIparams(1);
cy = ROIparams(2);
shapeFlag = ROIparams(3); % 0 = círculo, 1 = quadrado

hasRadius = numel(ROIparams) >= 4 && ~isempty(ROIparams(4));

if hasRadius
    r = ROIparams(4);
else
    r = [];
end

fillMode = isempty(nStims);
rMax = minDist * (1 + c);

% Gera pontos ao longo da tela inteira
xMin = 0; xMax = screenRes(1);
yMin = 0; yMax = screenRes(2);

% Funções de ROI
inCircle = @(x,y,rr) ((x-cx).^2 + (y-cy).^2) <= rr^2;
inSquare = @(x,y,rr) (abs(x-cx) <= rr) & (abs(y-cy) <= rr);

distToCircleEdge = @(x,y,rr) rr - sqrt((x-cx).^2 + (y-cy).^2);
distToSquareEdge = @(x,y,rr) rr - max(abs(x-cx), abs(y-cy));

% Grid para o Poisson disk sampling
cellSize = minDist / sqrt(2);
gridWidth  = ceil((xMax - xMin) / cellSize);
gridHeight = ceil((yMax - yMin) / cellSize);
accelGrid = cell(gridHeight, gridWidth);

toGrid = @(x,y) [floor((y-yMin)/cellSize)+1, floor((x-xMin)/cellSize)+1];

k = 30;

% Ponto inicial aleatório na tela
p0 = [rand*(xMax-xMin)+xMin, rand*(yMax-yMin)+yMin];

allPoints = p0;
activeList = p0;

gIdx = toGrid(p0(1), p0(2));
accelGrid{gIdx(1), gIdx(2)} = 1;

% POISSON SAMPLING
while ~isempty(activeList)

    idx = randi(size(activeList,1));
    basePoint = activeList(idx,:);
    found = false;

    for i = 1:k
        rr = minDist * (1 + c*rand);
        theta = 2*pi*rand;
        cand = basePoint + [rr*cos(theta), rr*sin(theta)];

        if cand(1) < xMin || cand(1) > xMax || cand(2) < yMin || cand(2) > yMax
            continue;
        end

        gIdx = toGrid(cand(1), cand(2));

        if gIdx(1)<1 || gIdx(1)>gridHeight || gIdx(2)<1 || gIdx(2)>gridWidth
            continue;
        end

        % collision
        rStart = max(1, gIdx(1)-2);
        rEnd   = min(gridHeight, gIdx(1)+2);
        cStart = max(1, gIdx(2)-2);
        cEnd   = min(gridWidth,  gIdx(2)+2);

        collision = false;
        for rIdx = rStart:rEnd
            for cIdx = cStart:cEnd
                pIdx = accelGrid{rIdx,cIdx};
                if isempty(pIdx), continue; end

                neighbors = allPoints(pIdx,:);
                distsSq = sum((neighbors - cand).^2, 2);
                if any(distsSq < minDist^2)
                    collision = true;
                    break;
                end
            end
            if collision, break; end
        end

        if ~collision
            allPoints = [allPoints; cand]; %#ok<*AGROW>
            activeList = [activeList; cand];
            accelGrid{gIdx(1), gIdx(2)} = size(allPoints,1);
            found = true;
            break;
        end
    end

    if ~found
        activeList(idx,:) = [];
    end
end

%% Seleção do ROI
% Caso 1: se não há raio nem número de estímulos, o ROI terá o tamanho
%         máximo delimitado pela tela e pelo seu centro
if ~hasRadius && fillMode
    r = min([cx, cy, screenRes(1)-cx, screenRes(2)-cy]);
end

% Caso 2: se não há raio mas há número de estímulos, expande ROI até 
%         conter no mínimo nStims estímulos
if ~hasRadius && ~fillMode
    r = minDist;
    step = minDist;

    while true
        if shapeFlag == 0
            inside = inCircle(allPoints(:,1), allPoints(:,2), r);
        else
            inside = inSquare(allPoints(:,1), allPoints(:,2), r);
        end

        if sum(inside) >= nStims
            break;
        end

        r = r + step;

        % segurança: não ultrapassar tela
        if r > max(screenRes)
            break;
        end
    end
end

% Caso 3: se há raio, seleciona apeas os pontos que etão dentro dele, e a 
%         uma distância mínima de suas bordas (para que o estímulo inteiro,
%         de raio rStim, caiba no ROI)
if hasRadius
    if shapeFlag == 0
        inside = inCircle(allPoints(:,1), allPoints(:,2), r);
        distEdge = distToCircleEdge(allPoints(:,1), allPoints(:,2), r);
    else
        inside = inSquare(allPoints(:,1), allPoints(:,2), r);
        distEdge = distToSquareEdge(allPoints(:,1), allPoints(:,2), r);
    end

    valid = inside & (distEdge >= rStim);
    
    roiPoints = allPoints(valid,:);
else
    roiPoints = allPoints;
end


%% Seleção dos pontos
% Caso 1: se há restrição no número de pontos, seleciona os nStims pontos 
%         no ROI mais próximos ao centro
if ~fillMode
    
    if hasRadius && size(roiPoints,1) < nStims
        error('Não foi possível encontrar pontos suficientes dentro do ROI.');
    end

    % distância ao centro do ROI
    dists = sqrt((roiPoints(:,1) - cx).^2 + (roiPoints(:,2) - cy).^2);
    
    % ordenar do mais interno → mais externo
    [~, order] = sort(dists, 'ascend');
    
    % selecionar os nStims mais internos
    selected = roiPoints(order(1:nStims), :);

else

% Caso 2: se NÃO há restrição no número de pontos, pega todos o spontos
%         dentro do ROI
    selected = roiPoints;
end

% Saída
stimCenters = selected';

if isempty(selected)
    fixCenter = [];
else
    fixCenter = selected(randi(size(selected,1)), :)';
end

end
function [fixCenter, stimCenters] = getStimLocations2(ROIparams, nStims, minDist, ~, ~)
    % Poisson disk sampling exato
    
    cx = ROIparams(1); cy = ROIparams(2);
    rx = ROIparams(3); ry = ROIparams(4);
    
    % Coordenadas para os 4 cantos do menor retângulo que contém a elipse
    xMin = cx - rx; xMax = cx + rx;
    yMin = cy - ry; yMax = cy + ry;
    
    % Função da elipse
    inEllipse = @(x,y) ((x-cx)/rx).^2 + ((y-cy)/ry).^2 <= 1;

    % Cria um grid oculto 
    cellSize = minDist / sqrt(2);
    gridWidth  = ceil((xMax - xMin) / cellSize);
    gridHeight = ceil((yMax - yMin) / cellSize);
    accelGrid = cell(gridHeight, gridWidth);
    
    toGrid = @(x,y) [floor((y-yMin)/cellSize)+1, floor((x-xMin)/cellSize)+1];
    
    k = 60; % attempts per active point
    
    % --- Step 1: initial point (fixation) ---
    while true
        p0 = [rand*(xMax-xMin)+xMin, rand*(yMax-yMin)+yMin];
        if inEllipse(p0(1), p0(2))
            break;
        end
    end
    
    allPoints = p0;
    activeList = p0;
    
    gIdx = toGrid(p0(1), p0(2));
    accelGrid{gIdx(1), gIdx(2)} = 1;
    
    % --- Step 2: grow samples ---
    while ~isempty(activeList) && size(allPoints,1) < (nStims + 1)
        
        idx = randi(size(activeList,1));
        basePoint = activeList(idx,:);
        
        found = false;
        
        for i = 1:k
            r = minDist * (1 + rand); % r = minDist * (1 + 0.05*rand);
            theta = 2*pi*rand;
            cand = basePoint + [r*cos(theta), r*sin(theta)];
            
            if ~inEllipse(cand(1), cand(2))
                continue;
            end
            
            gIdx = toGrid(cand(1), cand(2));
            
            if gIdx(1)<1 || gIdx(1)>gridHeight || gIdx(2)<1 || gIdx(2)>gridWidth
                continue;
            end
            
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
                allPoints = [allPoints; cand];
                activeList = [activeList; cand];
                
                accelGrid{gIdx(1), gIdx(2)} = size(allPoints,1);
                found = true;
                break;
            end
        end
        
        if ~found
            activeList(idx,:) = []; % remove exhausted point
        end
    end
    
    fixCenter = allPoints(1,:)';
    stimCenters = allPoints(2:end,:)';
end
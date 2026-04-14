function [fixCenter, stimCenters, rMax] = getStimLocations2(ROIparams, nStims, minDist, specialFix, skipEllipse, c)

    % Poisson disk sampling com busca restrita ao intervalo [minDist, minDist*(1+c)]
    % Para o algoritmo original, c = 1; para distancias mais uniformes,
    % usar c = 0.05

    if nargin < 5, skipEllipse = false; end
    if nargin < 6, c = 1; end

    % --- NEW: fill mode ---
    fillMode = isempty(nStims);

    rMax = minDist * (1 + c);
    
    cx = ROIparams(1); cy = ROIparams(2);
    rx = ROIparams(3); ry = ROIparams(4);
    
    xMin = cx - rx; xMax = cx + rx;
    yMin = cy - ry; yMax = cy + ry;
    
    inEllipse = @(x,y) ((x-cx)/rx).^2 + ((y-cy)/ry).^2 <= 1;
    inRect = @(x,y) (x >= xMin) & (x <= xMax) & (y >= yMin) & (y <= yMax);

    % --- Grid ---
    cellSize = minDist / sqrt(2);
    gridWidth  = ceil((xMax - xMin) / cellSize);
    gridHeight = ceil((yMax - yMin) / cellSize);
    accelGrid = cell(gridHeight, gridWidth);
    
    toGrid = @(x,y) [floor((y-yMin)/cellSize)+1, floor((x-xMin)/cellSize)+1];
    
    k = 60;
    
    % --- Step 1: initial point ---
    while true
        p0 = [rand*(xMax-xMin)+xMin, rand*(yMax-yMin)+yMin];
        if skipEllipse
            break
        else
            if inEllipse(p0(1), p0(2)), break; end
        end
    end
    
    allPoints = p0;
    activeList = p0;
    
    gIdx = toGrid(p0(1), p0(2));
    accelGrid{gIdx(1), gIdx(2)} = 1;
    
    % --- Step 2: growth ---
    if ~fillMode
        if specialFix
            nStimsAux = nStims;
        else
            nStimsAux = nStims + 1;
        end
    end

    while ~isempty(activeList) && (fillMode || size(allPoints,1) < nStimsAux)
        
        idx = randi(size(activeList,1));
        basePoint = activeList(idx,:);
        
        found = false;
        
        for i = 1:k
            r = minDist * (1 + c*rand);
            theta = 2*pi*rand;
            cand = basePoint + [r*cos(theta), r*sin(theta)];
            
            % --- ROI check ---
            if skipEllipse
                if ~inRect(cand(1), cand(2)), continue; end
            else
                if ~inEllipse(cand(1), cand(2)), continue; end
            end
            
            gIdx = toGrid(cand(1), cand(2));
            
            if gIdx(1)<1 || gIdx(1)>gridHeight || gIdx(2)<1 || gIdx(2)>gridWidth
                continue;
            end
            
            % --- collision ---
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
            
            % --- accept ---
            if ~collision
                allPoints = [allPoints; cand];
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

    % --- Output ---
    if fillMode
        % no special fixation logic in fill mode
        fixCenter = [];
        stimCenters = allPoints';
        return;
    end

    if specialFix
        % (unchanged)
        fx = ((allPoints(:,1)-cx)/rx).^2 + ((allPoints(:,2)-cy)/ry).^2;
        distToBoundary = abs(1 - fx);
        
        [aux, order] = sort(distToBoundary);
        order(aux > .5) = [];
        candidates = order(1:min(6, numel(order)));
        
        i1 = candidates(randi(numel(candidates)));
        p1 = allPoints(i1,:);
        
        others = candidates(candidates ~= i1);
        
        if isempty(others)
            p2 = p1;
        else
            dists = sum((allPoints(others,:) - p1).^2, 2);
            [~, idxMin] = min(dists);
            p2 = allPoints(others(idxMin),:);
        end
        
        pmid = (p1 + p2)/2;
        
        dir = [(pmid(1)-cx)/rx^2, (pmid(2)-cy)/ry^2];
        dir(2) = 0;
        
        if norm(dir) < 1e-6
            dir = [sign(pmid(1)-cx), 0];
        else
            dir = dir / norm(dir);
        end
        
        r = minDist * (1 + c*rand);
        fixCenter = pmid + r * dir;
        
        stimCenters = allPoints';
        
    else
        fixCenter = allPoints(1,:)';
        stimCenters = allPoints(2:end,:)';
    end
end
function [fixCenter, stimCenters, rMax] = getStimLocations5(ROIparams, nStims, minDist, ~)
% Essa versão produz um cordão de Poisson ao sortear apenas da vzinhança
% daqueles que têm menos de 2 vizinhos...

    rMax = minDist;
    k = 30;
    
    cx = ROIparams(1); cy = ROIparams(2);
    rx = ROIparams(3); ry = ROIparams(4);
    
    xMin = cx - rx; xMax = cx + rx;
    yMin = cy - ry; yMax = cy + ry;
    
    inEllipse = @(x,y) ((x-cx)/rx).^2 + ((y-cy)/ry).^2 <= 1;

    % --- Grid ---
    cellSize = minDist / sqrt(2);
    gridW = ceil((xMax - xMin)/cellSize);
    gridH = ceil((yMax - yMin)/cellSize);
    grid = cell(gridH, gridW);
    
    toGrid = @(x,y) [floor((y-yMin)/cellSize)+1, floor((x-xMin)/cellSize)+1];

    % --- Storage ---
    maxPts = nStims * 3;
    allPoints = zeros(maxPts,2);
    neighborCount = zeros(maxPts,1);
    nPoints = 0;

    % --- Helpers ---
    function collision = hasCollision(pt)
        g = toGrid(pt(1), pt(2));
        if g(1)<1 || g(1)>gridH || g(2)<1 || g(2)>gridW
            collision = true; return;
        end
        
        r0 = max(1, g(1)-2); r1 = min(gridH, g(1)+2);
        c0 = max(1, g(2)-2); c1 = min(gridW, g(2)+2);
        
        collision = false;
        for rr = r0:r1
            for cc = c0:c1
                idx = grid{rr,cc};
                if isempty(idx), continue; end
                d2 = sum((allPoints(idx,:) - pt).^2,2);
                if any(d2 < minDist^2)
                    collision = true;
                    return;
                end
            end
        end
    end

    function idxs = getNeighbors(pt)
        if nPoints == 0
            idxs = [];
            return;
        end
        d = sqrt(sum((allPoints(1:nPoints,:) - pt).^2,2));
        idxs = find(d >= minDist & d <= rMax);
    end

    % --- Step 1: seed first point ---
    while true
        p0 = [rand*(xMax-xMin)+xMin, rand*(yMax-yMin)+yMin];
        if inEllipse(p0(1), p0(2)), break; end
    end
    
    nPoints = 1;
    allPoints(1,:) = p0;
    
    g = toGrid(p0(1), p0(2));
    grid{g(1), g(2)} = 1;

    activeList = 1;

    % --- Step 2: force second point (initial neighbor) ---
    while true
        theta = 2*pi*rand;
        r = minDist;
        cand = p0 + [r*cos(theta), r*sin(theta)];
        
        if ~inEllipse(cand(1), cand(2)), continue; end
        if hasCollision(cand), continue; end
        
        nPoints = 2;
        allPoints(2,:) = cand;
        
        g = toGrid(cand(1), cand(2));
        grid{g(1), g(2)} = 2;
        
        neighborCount(1) = 1;
        neighborCount(2) = 1;
        
        activeList = [1 2];
        break;
    end

    % --- Step 3: Bridson with degree-based activity ---
    while ~isempty(activeList) && nPoints < nStims
        
        baseIdx = activeList(randi(numel(activeList)));
        base = allPoints(baseIdx,:);
        found = false;
        
        for i = 1:k
            
            theta = 2*pi*rand;
            r = minDist + (rMax - minDist)*rand;
            cand = base + [r*cos(theta), r*sin(theta)];
            
            if ~inEllipse(cand(1), cand(2)), continue; end
            if hasCollision(cand), continue; end
            
            neigh = getNeighbors(cand);
            if isempty(neigh), continue; end
            
            % --- accept ---
            nPoints = nPoints + 1;
            allPoints(nPoints,:) = cand;
            
            g = toGrid(cand(1), cand(2));
            grid{g(1), g(2)} = nPoints;
            
            % --- update graph ---
            neighborCount(nPoints) = numel(neigh);
            neighborCount(neigh) = neighborCount(neigh) + 1;
            
            % --- update active list ---
            activeList(end+1) = nPoints;
            
            % remove satisfied nodes (degree >= 2)
            activeList = activeList(neighborCount(activeList) < 2);
            
            found = true;
            break;
        end
        
        if ~found
            activeList(activeList == baseIdx) = [];
        end
    end

    % --- Trim ---
    allPoints = allPoints(1:nPoints,:);
    neighborCount = neighborCount(1:nPoints);
    
    if nPoints > nStims
        [~, order] = sort(neighborCount, 'descend');
        allPoints = allPoints(order(1:nStims),:);
    end

    stimCenters = allPoints';

    % --- Step 4: fixation ---
    centroid = mean(allPoints,1);
    
    while true
        theta = 2*pi*rand;
        r = minDist * (1 + 0.2*rand);
        cand = centroid + [r*cos(theta), r*sin(theta)];
        
        if ~inEllipse(cand(1), cand(2)), continue; end
        
        d = sqrt(sum((allPoints - cand).^2,2));
        if any(abs(d - minDist) < 0.25*minDist)
            fixCenter = cand';
            break;
        end
    end

end
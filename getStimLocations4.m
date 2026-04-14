function [fixCenter, stimCenters, rMax] = getStimLocations4(ROIparams, nStims, minDist, ~)

    % --- Parameters ---
    c = 0.25;                 % thickness of distance band
    rMax = minDist * (1 + c);
    k = 80;                   % attempts per iteration
    maxGlobalAttempts = 5000; % safety
    
    cx = ROIparams(1); cy = ROIparams(2);
    rx = ROIparams(3); ry = ROIparams(4);
    
    xMin = cx - rx; xMax = cx + rx;
    yMin = cy - ry; yMax = cy + ry;
    
    inEllipse = @(x,y) ((x-cx)/rx).^2 + ((y-cy)/ry).^2 <= 1;

    % --- Grid acceleration ---
    cellSize = minDist / sqrt(2);
    gridW = ceil((xMax - xMin) / cellSize);
    gridH = ceil((yMax - yMin) / cellSize);
    grid = cell(gridH, gridW);
    
    toGrid = @(x,y) [floor((y-yMin)/cellSize)+1, floor((x-xMin)/cellSize)+1];

    % --- Storage ---
    allPoints = zeros(nStims,2);
    neighborCount = zeros(nStims,1);
    nPoints = 0;

    % --- Helper: get neighbors ---
    function idxs = getNeighbors(pt)
        if nPoints == 0
            idxs = [];
            return;
        end
        d = sqrt(sum((allPoints(1:nPoints,:) - pt).^2,2));
        idxs = find(d >= minDist & d <= rMax);
    end

    % --- Helper: collision check ---
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

    % --- Step 1: initial seed ---
    while true
        p = [rand*(xMax-xMin)+xMin, rand*(yMax-yMin)+yMin];
        if inEllipse(p(1), p(2))
            nPoints = 1;
            allPoints(1,:) = p;
            g = toGrid(p(1), p(2));
            grid{g(1), g(2)} = 1;
            break;
        end
    end

    % --- Step 2: grow points ---
    attempts = 0;
    
    while nPoints < nStims && attempts < maxGlobalAttempts
        attempts = attempts + 1;
        
        % --- prioritize weak nodes ---
        weak = find(neighborCount(1:nPoints) < 2);
        if ~isempty(weak)
            baseIdx = weak(randi(numel(weak)));
        else
            baseIdx = randi(nPoints);
        end
        
        base = allPoints(baseIdx,:);
        accepted = false;
        
        for i = 1:k
            
            theta = 2*pi*rand;
            r = minDist * (1 + c*rand);
            cand = base + [r*cos(theta), r*sin(theta)];
            
            if ~inEllipse(cand(1), cand(2)), continue; end
            if hasCollision(cand), continue; end
            
            neigh = getNeighbors(cand);
            
            if isempty(neigh)
                continue;
            end
            
            % --- enforce ≥2 neighbors ---
            if numel(neigh) < 2
                if ~any(neighborCount(neigh) < 2)
                    continue;
                end
            end
            
            % --- accept ---
            nPoints = nPoints + 1;
            allPoints(nPoints,:) = cand;
            
            % update grid
            g = toGrid(cand(1), cand(2));
            grid{g(1), g(2)} = nPoints;
            
            % update graph
            neighborCount(nPoints) = numel(neigh);
            neighborCount(neigh) = neighborCount(neigh) + 1;
            
            accepted = true;
            break;
        end
        
        % --- small relaxation if stuck ---
        if ~accepted && attempts > maxGlobalAttempts/2
            c = min(c + 0.05, 0.6);
            rMax = minDist * (1 + c);
        end
    end

    % --- Step 3: repair (ensure all have ≥2 neighbors) ---
    for iter = 1:200
        weak = find(neighborCount(1:nPoints) < 2);
        if isempty(weak), break; end
        
        idx = weak(randi(numel(weak)));
        base = allPoints(idx,:);
        
        for i = 1:30
            theta = 2*pi*rand;
            r = minDist * (1 + 0.3*rand);
            cand = base + [r*cos(theta), r*sin(theta)];
            
            if ~inEllipse(cand(1), cand(2)), continue; end
            if hasCollision(cand), continue; end
            
            neigh = getNeighbors(cand);
            if isempty(neigh), continue; end
            
            nPoints = nPoints + 1;
            if nPoints > size(allPoints,1)
                allPoints = [allPoints; zeros(10,2)];
                neighborCount = [neighborCount; zeros(10,1)];
            end
            
            allPoints(nPoints,:) = cand;
            g = toGrid(cand(1), cand(2));
            grid{g(1), g(2)} = nPoints;
            
            neighborCount(nPoints) = numel(neigh);
            neighborCount(neigh) = neighborCount(neigh) + 1;
            break;
        end
    end

    % --- Trim to exactly nStims (keep well-connected ones) ---
    allPoints = allPoints(1:nPoints,:);
    neighborCount = neighborCount(1:nPoints);
    
    if nPoints > nStims
        [~, order] = sort(neighborCount, 'descend');
        allPoints = allPoints(order(1:nStims),:);
    end

    stimCenters = allPoints';

    % --- Step 4: fixation (independent) ---
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
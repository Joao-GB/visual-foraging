function points = getHexDistribution(center, dist, nRings, rotationDeg)
    % center: [x, y] coordinates of the central point
    % dist:   Distance between adjacent points (side length of hexagon)
    % nRings: Number of layers around center (0=1 point, 1=7 points, 2=19 points...)
    % rotationDeg: Rotation of the entire grid in degrees (counter-clockwise)

    % 1. Generate the Grid in "Axial Coordinates" (q, r)
    % We use a standard axial coordinate system for hex grids.
    % q = column, r = row (diagonal). 
    % Convert to Cartesian: x = dist * (3/2 * q), y = dist * (sqrt(3)/2 * q + sqrt(3) * r)
    % BUT, for simple ring generation, it's easier to just use basis vectors.
    if size(center,1) > size(center,2), center = center'; end
    
    % Basis vectors for moving along the hex grid axes (0, 60, 120, 180, 240, 300 deg)
    % We start at 0 degrees and step by 60.
    thetas = (0:60:300)' * (pi/180); 
    basisVecs = [cos(thetas), sin(thetas)] * dist;
    
    % Pre-allocate size: 1 + sum(6*i for i=1:nRings)
    totalPoints = 1 + 3 * nRings * (nRings + 1);
    relPoints = zeros(totalPoints, 2);
    
    % Ring 0 is just [0,0]
    count = 1;
    relPoints(1, :) = [0, 0];
    
    % Generate Rings 1 to nRings
    for r = 1:nRings
        % Start the ring: Move 'r' steps along the 240-degree vector (index 5)
        % to get to the "bottom left" corner of the ring, then walk around.
        % (Arbitrary start, but consistent). 
        % A simpler way: Start at vector 4 (180 deg) * r, then walk around.
        
        % Let's use the standard "Spiral" algorithm:
        % 1. Start at 'r' steps in direction 4 (West/180 deg) -> actually let's use direction 5 (240 deg) 
        %    to match standard "pointy top" or "flat top" conventions, but since we rotate later,
        %    it doesn't matter. Let's start "West" (-dist*r, 0).
        
        currPos = basisVecs(4,:) * r; % Move r steps West
        
        % Walk along the 6 directions
        % Order of directions to walk a full circle: 
        % If we started West, we go: NE, E, SE, SW, W, NW (approx)
        % The standard basis order is 0, 60, 120, 180, 240, 300.
        % If we start at 4 (180), the walk sequence is: 5, 0, 1, 2, 3, 4 ?? 
        % Actually, the sequence to trace a ring is:
        % Move direction 6 (300 deg), then 1 (0 deg), then 2 (60 deg)...
        
        % Let's do the "Vector Walk" explicitly:
        % Start at: vector(5) * r [which is 240 deg]
        currPos = basisVecs(5,:) * r; 
        
        % For each of the 6 sides of the hexagon ring...
        for side = 1:6
            % In this side, we take 'r' steps
            % The direction to walk depends on which side we are on.
            % Side 1 walks dir 1 (0 deg), Side 2 walks dir 2 (60 deg)...
            walkDir = basisVecs(side, :); 
            
            for step = 1:r
                currPos = currPos + walkDir;
                count = count + 1;
                relPoints(count, :) = currPos;
            end
        end
    end
    
    % 2. Rotation (Rigid Body)
    if rotationDeg ~= 0
        theta = rotationDeg * (pi/180);
        R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        % Apply rotation: (R * p')'
        relPoints = (R * relPoints')';
    end
    
    % 3. Translation
    points = relPoints + center;
end
function [nbhd1, nbhd2, nbhdElse] = getHexNeighborhoods(points, Pidx)
    if size(points,2) > size(points,1), points = points'; end

    if isempty(Pidx)
        aux = 1:size(points,2);
        nbhd1 = aux; nbhd2 = aux; nbhdElse = aux;
        return
    end
    
    % points: Nx2 matrix of all grid points
    % Pidx:   Scalar index of the center point P in the 'points' list
    
    % 1. Get coordinates of P
    P = points(Pidx,:);
    
    % 2. Calculate Euclidean distances to all other points
    % Vectorized calculation for speed
    dists = sqrt(sum((points - P).^2, 2));
    
    % 3. Determine the Grid Spacing (d)
    % We sort distances to find the closest non-zero distance.
    sortedDists = sort(dists);
    
    % Sanity check: Ensure we have at least 2 points
    if length(sortedDists) < 2
        nbhd1 = []; nbhd2 = []; nbhdElse = [];
        return;
    end
    
    % The first distance is 0 (point to itself). The second is 'd'.
    gridSpacing = sortedDists(2);
    
    % 4. Classify Neighbors by Distance Bands
    % We use thresholds halfway between the known geometric layers to be robust.
    % Ring 1 is at 1.0d.
    % Ring 2 is between 1.73d and 2.0d.
    % Ring 3 starts at 2.64d.
    
    % Threshold 1: Between Ring 1 (1.0) and Ring 2 (1.73). Let's use 1.4.
    limit1 = 1.4 * gridSpacing;
    
    % Threshold 2: Between Ring 2 (2.0) and Ring 3 (2.64). Let's use 2.3.
    limit2 = 2.3 * gridSpacing;
    
    % --- Get Indices ---
    % Nbhd1: Neighbors strictly within the first limit (excluding P itself)
    nbhd1 = find(dists > 0.1*gridSpacing & dists < limit1)';
    
    % Nbhd2: Neighbors between Ring 1 and Ring 3
    nbhd2 = find(dists > limit1 & dists < limit2)';
    
    % NbhdElse: Everything further away
    nbhdElse = find(dists > limit2)';
end
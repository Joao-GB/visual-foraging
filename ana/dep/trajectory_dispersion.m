function [diameter, farthest_pair] = trajectory_dispersion(x, y)
% TRAJECTORY_DISPERSION - Computes the dispersion (maximal pairwise distance) of a 2D trajectory.
%
% Inputs:
%   x : Vector of x-coordinates.
%   y : Vector of y-coordinates (same length as x).
%
% Outputs:
%   diameter       : Maximal Euclidean distance between any two points.
%   farthest_pair  : Indices of the two points defining this distance.

    if nargin < 2
        y = x(2,:);
        x = x(1,:);
    end
    
    % Validate inputs
    if length(x) ~= length(y)
        error('x and y must have the same length.');
    end
    
    % Combine coordinates into Nx2 matrix
    points = [x(:), y(:)];
    n = size(points, 1);
    
    % Initialize
    max_distance = 0;
    idx1 = 1;
    idx2 = 1;
    
    % Brute-force search for farthest pair (O(n^2))
    for i = 1:n
        for j = i+1:n
            dist = norm(points(i,:) - points(j,:));
            if dist > max_distance
                max_distance = dist;
                idx1 = i;
                idx2 = j;
            end
        end
    end
    
    % Assign outputs
    diameter = max_distance;
    farthest_pair = [idx1, idx2];
end
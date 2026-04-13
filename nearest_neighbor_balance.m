function metrics = nearest_neighbor_balance(points, tol)
% points: N x 2 matrix of [x y] coordinates
% tol: relative-difference threshold, e.g. 0.10 = 10%
%
% Returns a table with:
%   d1, d2            = nearest and second-nearest distances
%   var_sample        = sample variance of [d1 d2]
%   var_population    = population variance of [d1 d2]
%   abs_diff          = |d2 - d1|
%   rel_diff          = |d2 - d1| / mean(d1,d2)
%   pass              = true if rel_diff <= tol

    if nargin < 2
        tol = 0.10; % default: 10%
    end

    N = size(points, 1);

    % Pairwise Euclidean distance matrix
    dx = points(:,1) - points(:,1).';
    dy = points(:,2) - points(:,2).';
    D = sqrt(dx.^2 + dy.^2);

    % Ignore self-distances
    D(1:N+1:end) = Inf;

    % Sort each row: first two entries are nearest neighbors
    Ds = sort(D, 2, 'ascend');
    d1 = Ds(:,1);
    d2 = Ds(:,2);

    % Variance of the two nearest-neighbor distances
    var_sample     = var([d1 d2], 0, 2); % MATLAB default: divide by n-1
    var_population = var([d1 d2], 1, 2); % divide by n

    % More interpretable metrics
    abs_diff = abs(d2 - d1);
    rel_diff = abs_diff ./ ((d1 + d2)/2);

    % Criterion
    pass = rel_diff <= tol;

    metrics = table((1:N).', d1, d2, var_sample, var_population, ...
                    abs_diff, rel_diff, pass, ...
        'VariableNames', {'point','d1','d2','var_sample','var_population', ...
                          'abs_diff','rel_diff','pass'});
end
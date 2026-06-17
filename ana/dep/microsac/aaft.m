function xs = aaft(x);
%
%
%
%
N = length(x);
% 1. Ranking
[Sx, Rx] = sort(x);
% 2. Random Gaussian data
g = randn(1,N);
% 3. Sort Gaussian data
Sg = sort(g);
% 4. Rearrange Gaussian data
y(Rx) = Sg;
% 5. Create phase-randomized surrogate
y1 = ftpr(y);
% 6. Ranked time series
[Sy, Ry1] = sort(y1);
% 7. AAFT surrogate time series
xs(Ry1) = Sx;

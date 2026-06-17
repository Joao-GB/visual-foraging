function y = rmvsac(x,sac);
%
%
%

% 1. List saccade samples
del = [];
for i=1:size(sac,1)
    del = [del sac(i,1):sac(i,2)];
end
del = del';

% 2. Remove velocity samples
x0 = x(1,:);
dx = diff(x);
dx(del,:) = [];
dx = [x0; dx];

% 3. Cumulative sm
y = cumsum(dx);


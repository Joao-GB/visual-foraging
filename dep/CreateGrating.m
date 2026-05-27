function [patch] = CreateGrating(size,frequency,angle,phase,contrast)

if nargin < 5
    error('Wrong input argument list.');
end

size = floor(size/2)*2;
[x,y] = meshgrid((1:size)-(size+1)/2);

patch = 0.5*contrast*cos(2*pi*(frequency*(sin(pi/180*angle)*x+cos(pi/180*angle)*y)+phase));
end
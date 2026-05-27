function [a, b]= dva2pix(dist,width,res,dva)
% [pixels of dva, pixelsize]= dva2pix(dist,width,res,dva)
%
% calculates de number of pixels of a given degree of visual angle

% dist = distance to the screen in cm
% width = width of the screen in cm
% res = screen resolution in pixels (e.g. 1080)
% dva = degrees of visual angle diameter

b = width/res; %pixel size in cm 

sz = 2*dist*tan(pi*dva/360);
a = round(sz/b);
end


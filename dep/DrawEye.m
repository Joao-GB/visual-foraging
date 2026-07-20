function DrawEye(win, center, radius, irisColor, alphaColor)
% Desenha um olho na janela win, com centro, raio e cor da írisespecificados
if nargin < 5, alphaColor = .72; end
x = center(1);
y = center(2);

%% Esclera (com um cortorno escuro)
rx = radius * 2.;
ry = radius * 1.3;

irisColor = [irisColor alphaColor];

t = linspace(0, pi, 40);

upper = [x + rx*cos(t)', ...
         y - ry*(sin(t).^1.2)'];

lower = [x + rx*cos(fliplr(t))', ...
         y + ry*(sin(fliplr(t)).^1.2)'];

pts = [upper; lower];

Screen('FillPoly', win, [1 1 1, alphaColor], pts);
Screen('FramePoly', win, [0 0 0, alphaColor], pts, 3);

%% Íris
irisRect = CenterRectOnPoint([0 0 2.2*radius 2.2*radius], x, y);
Screen('FillOval', win, irisColor, irisRect);

%% Pupila
pupilRadius = radius * 0.45;
pupilRect = CenterRectOnPoint([0 0 2*pupilRadius 2*pupilRadius], x, y);
Screen('FillOval', win, [0 0 0 alphaColor], pupilRect);

end
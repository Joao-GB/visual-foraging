function DrawCalibrationTarget(win, x, y, varargin)
% rect = DrawCalibrationTarget(win, x, y)
% rect = DrawCalibrationTarget(win, x, y, el)
% rect = DrawCalibrationTarget(win, x, y, sizePct, widthPct, fg, bg)
%
%
% Inputs
%   win      PTB window pointer
%   x,y      Target center
%
% Optional
%   el       EyeLink defaults structure
%
% or
%
%   sizePct  Outer target diameter as % of screen width (default 2.5)
%   widthPct Width of the central hole as % of screen width (default 0.7)
%   fg       Outer color (default white)
%   bg       Inner color (default background/black)
%
% Output
%   rect     Bounding rectangle of the outer target

% Defaults
sizePct  = 2;
widthPct = 0.3;
fg = 0;
bg = .5;

if nargin >= 4
    if isstruct(varargin{1})
        el = varargin{1};
        sizePct  = el.calibrationtargetsize;
        widthPct = el.calibrationtargetwidth;
        fg = el.calibrationtargetcolour;
        bg = el.backgroundcolour;
    else
        if numel(varargin) >= 1, sizePct  = varargin{1}; end
        if numel(varargin) >= 2, widthPct = varargin{2}; end
        if numel(varargin) >= 3, fg = varargin{3}; end
        if numel(varargin) >= 4, bg = varargin{4}; end
    end
end

[w, ~]=Screen('WindowSize', win);

outer = round(sizePct  /100 * w);
inner = round(widthPct /100 * w);

Screen('FillOval', win, fg, ...
    [x-outer/2 y-outer/2 x+outer/2 y+outer/2]);

Screen('FillOval', win, bg, ...
    [x-inner/2 y-inner/2 x+inner/2 y+inner/2]);
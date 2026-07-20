function DrawSlideExample(win, rect, leftTex, rightTex)

% Screen dimensions
[w, h] = RectSize(rect);

%% Layout parameters
margin      = 40;
titleY      = 40;

imageWidth  = round(0.30*w);
imageHeight = round(0.40*h);

imageY = round(0.18*h);

leftX  = round(w*0.25);
rightX = round(w*0.75);

%% Title
Screen('TextSize', win, 28);
DrawFormattedText(win, ...
    'Example stimuli', ...
    margin, titleY, 1);

%% Destination rectangles
leftRect = CenterRectOnPointd(...
    [0 0 imageWidth imageHeight], ...
    leftX, imageY + imageHeight/2);

rightRect = CenterRectOnPointd(...
    [0 0 imageWidth imageHeight], ...
    rightX, imageY + imageHeight/2);

%% Draw textures
Screen('DrawTexture', win, leftTex, [], leftRect);
Screen('DrawTexture', win, rightTex, [], rightRect);

%% Labels
Screen('TextSize', win, 20);
Screen('TextSize', dpP.window, prm.textSizeNormal);
DrawFormattedText(win, ...
    'Look to the LEFT image.', ...
    leftRect(1), leftRect(4)+20, 1);

DrawFormattedText(win, ...
    'Look to the RIGHT image.', ...
    rightRect(1), rightRect(4)+20, 1);

end
function foragingFlip1(win, centers, dstCoord, idxForOvals, s, colors, aT, tOri, pW, fb, wrongCol, rightCol)
    if nargin >= 10
        rightCount = sum(fb == 1);
        wrongCount = sum(fb == 0);
        colors(:, fb == 1) = repmat(rightCol', [1, rightCount]);
        colors(:, fb == 0) = repmat(wrongCol', [1, wrongCount]);
    end
    Screen('FrameOval', win, colors(:, idxForOvals), dstCoord(:, idxForOvals), pW(:, idxForOvals));

    targets = (aT == 1);
%     neutralStim = (aT == -1);
    nonTargets = (aT == 0);
    drawInclinedLines(win, centers(:,targets), s*.8, tOri, colors(:,targets), pW(:, targets));
    drawDots(win, centers(:,nonTargets), s/5, colors(:,nonTargets));
end


function drawDots(windowPtr, centers, len, color)
    if ~isempty(centers)
        Screen('DrawDots', windowPtr, centers, len, color, [], 2);
    end
end


function drawInclinedLines(windowPtr, centers, len, ang, color, width)

    if ~isempty(centers)
        % Simplified version for uniform lines
        numLines = size(centers, 2);
        angRad = deg2rad(ang);
        hl = len/2;
        
        % Calculate offsets
        dx = hl .* cos(angRad);
        dy = hl .* sin(angRad);
        
        % Create interleaved start and end points
        xy = [centers(1,:) - dx; centers(2,:) - dy; 
              centers(1,:) + dx; centers(2,:) + dy];
        xy = reshape(xy, 2, 2 * numLines);
        
        Screen('DrawLines', windowPtr, xy, width, color(:,repelem(1:size(color,2), 2)));
    end
end
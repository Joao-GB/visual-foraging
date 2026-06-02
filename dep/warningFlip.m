function currTime = warningFlip(win, centers, dstCoord, idx, s, pW, warningCol)
    colors = repmat(warningCol', [1, numel(idx)]);
    Screen('FrameOval', win, colors, dstCoord(:, idx), pW(:, idx));

    drawDots(win, centers(:,idx), s/2, colors, 1);

    currTime = Screen('Flip', win);
    WaitSecs(0.001); 
end


function drawDots(windowPtr, centers, len, color, contour)
if nargin < 5, contour = false; end
    if ~isempty(centers)
        if contour
            Screen('DrawDots', windowPtr, centers, len, [0 0 0], [], 2);
        end
        Screen('DrawDots', windowPtr, centers, len*(3/4), color, [], 2);
    end
end
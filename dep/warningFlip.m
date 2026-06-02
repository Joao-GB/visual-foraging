function currTime = warningFlip(win, centers, dstCoord, idx, s, pW, warningCol)
    colors = repmat(warningCol', [1, numel(idx)]);
    Screen('FrameOval', win, colors(:, idx), dstCoord(:, idx), pW(:, idx));

    drawDots(win, centers(:,idx), s/5, colors(:, idx));

    currTime = Screen('Flip', win);
    WaitSecs(0.001); 
end


function drawDots(windowPtr, centers, len, color)
    if ~isempty(centers)
        Screen('DrawDots', windowPtr, centers, len, color, [], 2);
    end
end
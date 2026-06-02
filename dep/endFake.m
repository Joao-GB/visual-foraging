function endFake(aux, dpP, drP, prm)
%% Faz uma transição lenta à cor original
    if aux.logoTex ~= -1, Screen('Close', aux.logoTex); end
    Screen('Close', aux.blockTex);

    img = Screen('GetImage', dpP.window);
    frozenTex = Screen('MakeTexture', dpP.window, img);

    fadeStart = GetSecs;
    while true
        elapsed = GetSecs - fadeStart;
        t = min(elapsed / prm.fadeInDur1, 1);
        texAlpha = 1-t;
    
        Screen('FillRect', dpP.window, drP.grey);
        Screen('DrawTexture', dpP.window, frozenTex, [], [], [], [], texAlpha);
    
        Screen('Flip', dpP.window);
    
        if t >= 1
            break;
        end
    end
    
    Screen('Close', frozenTex);
    WaitSecs(prm.fadeInDelay1);
end
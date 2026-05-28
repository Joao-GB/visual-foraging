function decision = pauseScreen(tkP, dpP, drP, prm, pauseMode)
    decision = '';

    leaveOptions = {'quit', 'gotoMenu'};
    parentDir = prm.currFolder;
    leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3};
    cX = dpP.winCenter(1); cY = dpP.winCenter(2);


    titleText = 'PAUSADO'; titleMargin = cX*prm.titleMarginFactor;
    %% Tipos de botões
    switch pauseMode
        case 'menu'
            options     = {'quit', 'resume'};
            optionsName = {'Sair', 'Prosseguir'};

        case 'trial'
            options = {'gotoMenu', 'restartBlock', 'recalibrate', 'resume'};
            optionsName = {'Menu', 'Reiniciar Bloco', 'Recalibrar', 'Prosseguir'};

        case 'block'
            options = {'gotoMenu', 'recalibrate', 'resume'};
            optionsName = {'Menu', 'Recalibrar', 'Prosseguir'};
    end

    icons    = cellfun(@(x) fullfile(parentDir, prm.imgFolder, [x prm.imgExtension]), options, 'UniformOutput', false);
    iconsImg = cellfun(@(x) alphaRead(x), icons, 'UniformOutput', false);
    imgSizes = cellfun(@(x) size(x), iconsImg, 'UniformOutput', false);
    auxImg   = cellfun(@(x) drP.black*ones(x(1), x(2), 4), imgSizes, 'UniformOutput', false);
    iconsImg = cellfun(@(x, a) addAlpha(x,a), auxImg, iconsImg, 'UniformOutput', false);
    iconsTex = cellfun(@(x) Screen('MakeTexture', dpP.window, x), iconsImg);
    clear icons imgSizes auxImg iconsImg;

    L = numel(options);
    if (L+1)/2 == floor((L+1)/2)
        firstSelectL = (L+1)/2;
        firstSelectR = (L+1)/2;
    else
        firstSelectL = floor((L+1)/2);
        firstSelectR = ceil((L+1)/2);
    end
    selected = 0;

    %% Retângulos dos botões
    gap  = 80;
    btnW = min(300, (dpP.winRect(3)-(L+1)*gap)/L); btnH = 300;
    totalW = L*btnW + (L-1)*gap;
    leftX0 = cX - totalW/2;
    leftX  = leftX0 + (0:L-1) * (btnW + gap);

    btnRects = zeros(4, L);
    btnRects(1,:) = leftX; btnRects(2,:) = cY - btnH/2;
    btnRects(3,:) = leftX + btnW; btnRects(4,:) = cY + btnH/2;
    
    prm.upFrac = 2/3; SizeFactor = .6;
    RectSplit = btnRects(2,:)+btnH*prm.upFrac;
    upRects = btnRects; upRects(4,:) = RectSplit; downRects = btnRects; downRects(2,:) = RectSplit;
    uRH = RectHeight(upRects(:,1)'); uRW = RectWidth(upRects(:,1)');
    uRHW = min(uRH, uRW);
    for i=1:L
        [cXaux, cYaux] = RectCenter(upRects(:,i));
        upRects(:,i) = CenterRectOnPoint([0 0, uRHW, uRHW]*SizeFactor, cXaux, cYaux);
    end
    clear cXaux cYaux

    keyWasDown = false;
    lastActionTime = 0; lastRepeatTime = -Inf;
    doAction = false;

    KbReleaseWait;

    while true
        Screen('TextStyle', dpP.window, 0);
        Screen('TextSize', dpP.window, prm.textSizeTitle);
        Screen('DrawText', dpP.window, titleText, titleMargin, titleMargin,  drP.blackGrey);

        Screen('TextSize', dpP.window, prm.textSizeNormalish);

        %% Desenha os retângulos e textos e texturas
        Screen('TextStyle', dpP.window, 1);
        Screen('FillRect', dpP.window, drP.blue, btnRects);
        Screen('DrawTextures', dpP.window, iconsTex, [], upRects);

        for i=1:L
            DrawFormattedText(dpP.window, optionsName{i}, 'center', 'center',  drP.black, [], [], [], [], [], downRects(:,i)');
        end

        % Se houver algum retângulo selecionado, adiciona contorno
        if selected ~= 0
            Screen('FrameRect', dpP.window,  drP.darkBlue, btnRects(:, selected), prm.pW2);
        end
        Screen('Flip', dpP.window);

        [keyIsDown, ~, keyCode] = KbCheck;

        if keyIsDown
            if ~keyWasDown
                lastActionTime = GetSecs;
                doAction = true;
            elseif GetSecs - lastActionTime > prm.repeatDelay
                if GetSecs - lastRepeatTime > prm.repeatRate
                    doAction = true;
                    lastRepeatTime = GetSecs;
                end
            end
            if doAction
                Screen('TextStyle', dpP.window, 0);
                doAction = false;
                if keyCode(leftKey)
                    if selected == 0
                        selected = firstSelectL;
                    else
                        selected = max(1, selected - 1);
                    end
    
                elseif keyCode(rightKey)
                    if selected == 0
                        selected = firstSelectR;
                    else
                        selected = min(L, selected + 1);
                    end
    
                elseif keyCode(spaceKey)
                    if selected ~= 0
                        decision = options{selected};
                    end
                    KbReleaseWait
                    if ismember(decision, leaveOptions)
                        KbReleaseWait;
                        confirmed = confirmationScreen(dpP, drP, prm, tkP.keys);
                        if confirmed
                            return;
                        else
                            decision = '';
                        end
                    end
                    if ~isempty(decision)
                        for i=1:L, Screen('Close', iconsTex(i)); end
                        return;
                    end
                end
            end
        else
            lastRepeatTime = -Inf;
        end
        keyWasDown = keyIsDown;
        WaitSecs(0.01);
    end
end
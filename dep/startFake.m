function aux = startFake(tkP, dpP, drP, prm, loadingMode, txP, ori)

    leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; escapeKey = tkP.keys{4};

    parentDir = prm.currFolder;
    fileLogo1 = 'logo_bBg_wSb.png'; fileLogo2 = 'logo_gBg_wSb.png'; fileRoundRect = '1by2_rect.png';

    pathLogo1 = fullfile(parentDir, prm.imgFolder, fileLogo1);
    pathLogo2 = fullfile(parentDir, prm.imgFolder, fileLogo2);

    pathRectRound = fullfile(parentDir, prm.imgFolder, fileRoundRect);
    clear parentDir fileRoundRect

    if strcmp(loadingMode, 'opening')
        pathLogo = pathLogo1; f = 0;
        bgColor   = drP.black;
        totalTime = 15;
        hints     = {
                     'Ganhe brinde extra ao indicar outras pessoas\n para participarem de nossos experimentos!'
                     'Faça parte da lista de sujeitos experimentais\n do laboratório!'
                     'Evite distrações enquanto realiza a tarefa.'
                     'Se o experimento tiver mais de uma sessão, \nnão esqueça de comparecer em todas elas!'
                    };

        hintChangeInterval = 6; % Tempo durante o qual uma dica é exibida
        showBar     = true;
        barYOffset  = 135;
        hintYOffset = 75;       % Quantos pixels abaixo da barra aparecem as dicas
    elseif ismember(loadingMode, {'start', 'trial'})
        pathLogo = pathLogo2; f = 1;
        bgColor   = drP.grey;
        totalTime = 8;
        hints     = {
                     'Memorize apenas onde estão os alvos!'
                     'Faça pausas apenas quando o tipo de alvo mudar.'
                     'Evite se mexer durante o experimento.'
                     'Tente se manter concentrado e relaxado.'
                     'Priorize o desempenho, e não a velocidade.'
                     'Não esqueça de olhar para a cruz no início\nde cada tentativa!'
                    };
        hintChangeInterval = 4.5;
        showBar     = false;
        barYOffset  = 105;
        hintYOffset = 105;
    end
    
    %% Parâmetros
    nBlocks          = 20;      % Número de blocos na barra de progresso
    nSteps           = 20;      % São dados tantos passos quanto os blocos
    barHeight        = 20;      % Altura da barra de progresso em pixels
    blockGap         = 6;       % Espaço entre blocos em pixels
    imageYOffsetFrac = 0.15;    % Quão para cima ficará a imagem em relação ao centro
    creditXOffset    = 50;
    creditYOffset    = 25;
    arrowsXOffset    = 10;

    blockImg = imread(pathRectRound);
    
    white = repmat(drP.white, 1, 3);
    
    Screen('TextFont', dpP.window, prm.textFont);
    Screen('TextColor', dpP.window, drP.white);
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    cX = dpP.winCenter(1); cY = dpP.winCenter(2);
    
    %% Carrega a imagem do logo e ajusta seu tamanho (no máximo 2/3 da tela)
    
    logoImg = imread(pathLogo);
    % Cria canal alfa para remover o fundo (f é parâmetro ajustado conforme 
    % a cor do fundo)
    if size(logoImg,3) == 3
        meanLogo = mean(logoImg, 3);
        alpha = uint8( 255*(meanLogo-f*min(meanLogo(:)))/( max(meanLogo(:)) - min(meanLogo(:)) ) );
        logoImg(:,:,4) = alpha;
    end
    logoTex = Screen('MakeTexture', dpP.window, logoImg);
    clear pathLogo1 pathLogo2 pathLogo pathRectRound
    
    [imgH, imgW, ~] = size(logoImg);
    maxImgW = dpP.winRect(3) * 1/3; maxImgH = dpP.winRect(4) * 1/3;
    
    scaleFactor = min([1, maxImgW/imgW, maxImgH/imgH]);
    scaledW = imgW * scaleFactor; scaledH = imgH * scaleFactor;
    
    imgRect = CenterRectOnPointd([0 0 scaledW scaledH], cX, cY - dpP.winRect(4) * imageYOffsetFrac);
    clear logoImg imgH imgW maxImgW maxImgH scaleFactor scaledW scaledH

    
    if strcmp(loadingMode, 'trial')
        totalTime = totalTime*(3/8);
        Screen('Close', logoTex);
        logoTex = -1;
    end

    %% Retângulo dos créditos
    creditX = imgRect(1) - creditXOffset; creditY = imgRect(2) - creditYOffset;
    if logoTex == -1
        creditText = 'ORIENTAÇÃO DO ALVO:';
    else
        creditText = 'UM EXPERIMENTO DE';
    end
    
    %% Carrega a imagem do retângulo arredondado como textura
        % Se não tiver canal alfa, cria de modo que o preto seja
        % transparente (para usar a cor do fundo)
    if size(blockImg,3) == 3
        alpha = uint8(mean(blockImg, 3));
        blockImg(:,:,4) = alpha;
    end
    
    blockTex = Screen('MakeTexture', dpP.window, blockImg);
    
    %% Barra de progresso
    % Tamanho original da imagem do bloco
    [bH, bW, ~] = size(blockImg);
    blockAspect = bW / bH;
    
    % A altura da imagem é dada com base na altura da barra
    blockHeight = barHeight;
    blockWidth  = blockHeight * blockAspect;
    
    % A largura da barra é dada com base no tamanho dos blocos e no espaço
    % entre eles
    totalBlocksWidth = nBlocks * blockWidth + (nBlocks-1) * blockGap;
    barX = cX - totalBlocksWidth / 2;
    barY = cY + barYOffset;
    
    % Contorno da barra
    barOutline = [barX - blockGap, barY - blockGap, barX + totalBlocksWidth + blockGap, barY + barHeight + blockGap];
    [barCX, barCY] = RectCenter(barOutline);
    
    % Retângulos em que são apresentadas as imagens dos blocos
    blockRects = zeros(4, nBlocks);
    for i = 1:nBlocks
        x1 = barX + (i-1) * (blockWidth + blockGap);
        blockRects(:,i) = [x1 barY x1+blockWidth barY+blockHeight];
    end
    clear blockImg bH bW blockAspect blockHeight blockWidth totalBlocksWidth
    clear barHeight barX blockGap imageYOffsetFrac creditXOffset creditYOffset 
    
    % As durações de cada passo são sorteadas uniformemente entre 0 e 1, de
    % modo que a soma coincida com a duração total
    stepDurations = rand(1, nSteps);
    stepDurations = stepDurations / sum(stepDurations) * totalTime;
    stepTimes = cumsum(stepDurations);
    
    % Aumenta um pouquinho o tempo total para ter algum instante com o 10/10
    totalTime = totalTime + prm.repeatDelay;
    blocksPerStep = nBlocks / nSteps;
    
    %% Dicas
    Screen('TextSize', dpP.window, prm.textSizeNormal);
    hintY = barY + hintYOffset;
    maxW = 0; maxH = 0;
    for i = 1:numel(hints)
        [~, ~, hintsRect] = DrawFormattedText(dpP.window, hints{i});
        currW = RectWidth(hintsRect); currH = RectHeight(hintsRect);
        if currW > maxW, maxW = currW; end
        if currH > maxH, maxH = currH; end
    end
    hintsRect = CenterRectOnPointd([0, 0, maxW, maxH], cX, hintY);
    clear stepDurations maxW maxH currW currH

    %% Carregando alternativo
    nDots = 4; amp = 8; freq = 2; k = 3;
    phaseStep = 2*pi/nDots; phase = (0:(nDots-1))*phaseStep;
    adjSin = @(x) amp .* sin(2*pi*freq*x + phase) .* (mod(floor(freq*x), k) == 0) .*(0.5 - 0.5*cos(2*pi*mod(freq*x,1)));
    radius = 3; distance = 3;
    loadingText = 'CARREGANDO';
    Screen('TextSize', dpP.window, prm.textSizeBig);
    boundsLoading = Screen('TextBounds', dpP.window, loadingText);
    loadingRect = CenterRectOnPointd(boundsLoading, barCX, barCY);
    loadingWidth  = RectWidth(boundsLoading);
    boundsAux = Screen('TextBounds', dpP.window, 'A');
    Screen('TextSize', dpP.window, prm.textSizeNormal);

    dotsX = cX + loadingWidth/2 + cumsum(repmat(distance+radius, 1, nDots));
    dotsY = repmat(loadingRect(2) + RectHeight(boundsAux), 1, nDots);

    %% Setas
    lArrowText = '<'; rArrowText = '>';
    arrowsRect = Screen('TextBounds', dpP.window, lArrowText);
    arrowsWidth = RectWidth(arrowsRect);
    larrowRect = CenterRectOnPointd(arrowsRect, hintsRect(1)-arrowsWidth-arrowsXOffset, hintY);
    rarrowRect = CenterRectOnPointd(arrowsRect, hintsRect(3)+arrowsWidth+arrowsXOffset, hintY);
    clear boundsLoading loadingWidth boundsAux nDots arrowsXOffset arrowsRect
    
    %% Loop principal
    hintIdx = 1;
    hints = hints(randperm(numel(hints)));
    currentHint = hints{mod(hintIdx, numel(hints))+1};
    hintPhase = 'steady';   % 'fadein', 'steady', 'fadeout'
    steadyStart = 0;
    hintAlpha = 1;
    
    if showBar, filledBlocks = 0; else, filledBlocks = nBlocks; end
    startTime = GetSecs;
    stepIndex = 0;
    
    doAction = false; selected = 0; keyWasDown = 0;
    Screen('TextSize', dpP.window, prm.textSizeSmall);
    lastActionTime = GetSecs;
    lastRepeatTime = GetSecs;
    while true
        tNow = GetSecs;
        elapsedTime = tNow - startTime;
    
        % A tela é interrompida apenas se excedeu o limite de tempo e se
        % todos os blocos foram preenchidos
        if elapsedTime >= totalTime && filledBlocks == nBlocks
            break;
        end
    
        if showBar
            % Atualiza a quantidade de blocos preenchidos
            if stepIndex < nSteps && elapsedTime >= stepTimes(stepIndex + 1)
                stepIndex = stepIndex + 1;
            end
            filledBlocks = round(stepIndex * blocksPerStep);
        end

        %% Setas de seleção das dicas
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
            if doAction && strcmp(hintPhase, 'steady')
                doAction = false;
                if keyCode(leftKey)
                    selected = -1;
                    hintPhase = 'fadeout'; fadeStart = elapsedTime;
                elseif keyCode(rightKey)
                    selected = 1;
                    hintPhase = 'fadeout'; fadeStart = elapsedTime;
                elseif keyCode(escapeKey)
                    break;
                end
            end
        end
        keyWasDown = keyIsDown;
    
        %% Atualiza a dica exibida
        % Se a dica ficou por tempo suficiente, ela deve sair e vir outra
        switch hintPhase
            case 'fadein'
                t = (elapsedTime - fadeStart) / prm.fadeInDur1;
                hintAlpha = min(1, t);
                % Se acabou a fase de transição (normalizada), muda
                % para estado steady
                if t >= 1
                    hintPhase = 'steady';
                    steadyStart = elapsedTime;
                end
        
            case 'steady'
                hintAlpha = 1;
        
                if elapsedTime - steadyStart >= hintChangeInterval - 2*prm.fadeInDur1
                    hintPhase = 'fadeout';
                    fadeStart = elapsedTime;
                end
        
            case 'fadeout'
                t = (elapsedTime - fadeStart) / prm.fadeInDur1;
                hintAlpha = max(0, (1 - t));
        
                if t >= 1
                    hintPhase = 'fadein';
                    fadeStart = elapsedTime;
                    if selected ~= 0
                        hintIdx = hintIdx + selected;
                        selected = 0;
                    else
                        hintIdx = hintIdx + 1;
                    end
                    currentHint = hints{mod(hintIdx, numel(hints))+1};
                end
        end

        Screen('FillRect', dpP.window, bgColor);
    
        %% Desenha o logo e o texto logo acima
        
        if logoTex == -1
            Screen('TextSize', dpP.window, prm.textSizeNormal);
            Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
            Screen('DrawTexture', dpP.window, txP.exampleGabor.tex, [], [], ori, [], [], [], [], []);
            Screen('DrawTexture', dpP.window, txP.exampleNoise.tex, [], [], [], [], [], [], [], []);
            Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTextures', dpP.window, txP.exampleBlob.tex, [], [], ori, [], [], [0 0 0 1]', [], [], txP.exampleBlob.props); 
            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        else
            Screen('DrawTexture', dpP.window, logoTex, [], imgRect);
        end
        DrawFormattedText(dpP.window, creditText, creditX, creditY, drP.white);
    
        %% Desenha a barra de progresso
        Screen('TextSize', dpP.window, prm.textSizeNormal);
        if showBar
            Screen('FrameRect', dpP.window, drP.white, barOutline, 2);
            if filledBlocks > 0
                Screen('DrawTextures', dpP.window, blockTex, [], blockRects(:,1:filledBlocks));
            end
        
            % Texto das etapas concluídas
            progressText = sprintf('%d/10%', floor(10 * filledBlocks / nBlocks));
            DrawFormattedText(dpP.window, progressText, 'center', barY - 25, drP.white);
        else
            Screen('TextSize', dpP.window, prm.textSizeBig);
            DrawFormattedText(dpP.window, loadingText, 'center', 'center', drP.white, [], [], [], [], [], loadingRect);
            Screen('DrawDots', dpP.window, [dotsX; dotsY+adjSin(GetSecs - startTime)], radius, drP.white, []);
            Screen('TextSize', dpP.window, prm.textSizeNormal);
        end

        %% Desenha/escreve a dica
        DrawFormattedText(dpP.window, currentHint, 'center', 'center', [white hintAlpha], [], [], [], [], [], hintsRect);
%         Screen('FrameRect', dpP.window, [], hintsRect);

        %% Desenha as setas
        Screen('TextSize', dpP.window, prm.textSizeSmall);
        Screen('TextStyle', dpP.window, max(-1*selected, 0));
        DrawFormattedText(dpP.window, lArrowText, 'center', 'center', drP.white, [], [], [], [], [], larrowRect);
        Screen('TextStyle', dpP.window, max(selected, 0));
        DrawFormattedText(dpP.window, rArrowText, 'center', 'center', drP.white, [], [], [], [], [], rarrowRect);
        Screen('TextStyle', dpP.window, 0);
%         Screen('TextSize', dpP.window, prm.textSizeNormal);
    
        Screen('Flip', dpP.window);
    end
    aux.logoTex  = logoTex;
    aux.blockTex = blockTex;
end

function [tkP, taskState] = menuScreen1(tkP, dpP, drP, txP, debug, prm)
%% Interpretação da variável taskState:
    % Linhas: treino e experimento; 
    % Colunas: começou e concluiu
    taskState = [0 0;
                 0 0];

    Screen('TextFont', dpP.window, prm.textFont);
    Screen('TextSize', dpP.window, prm.textSizeNormal);
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Reposiciona o mouse
    cX = dpP.winCenter(1); cY = dpP.winCenter(2);
    SetMouse(cX, cY, dpP.window);
    if debug == 0,  HideCursor(dpP.window); end

    titleText = 'MENU'; titleMargin = cX*prm.titleMarginFactor;
    currentScreen = 'main'; % Pode ser 'main' ou 'training'
    refreshLayout = true;

    % Cria texturas de cursor e treino
    mainOptions     = {'staircase', 'training', 'experiment'};
    mainOptionsName = {'Staircase', 'Treino', 'Experimento'};
    mainOptionsMsg  = {['Ajusta os parâmetros do ruído conforme o\n'...
               'participante usando o método de staircase.'],...
               ['Diferentes tipos de treino para o experimento,\n'...
               'com e sem movimentos oculares, e feedback ao \n'...
               'fim de cada tentativa.\n\n'...
               'Recomendado antes de toda sessão experimental.'],...
               'Inicia a sessão experimental propriamente.'};

    trainOptions     = {'cursor', 't1', 't2'};
    trainOptionsName = {'Demo Mouse', 'Treino 1', 'Treino 2'};
    trainOptionsMsg  = {'1', '2', '3'};
    
    skipText = 'Voltar';
    skipRect = Screen('TextBounds', dpP.window, skipText); 
    skipRect(1) = 60;                                     % Posição X (Canto esquerdo)
    skipRect(2) = dpP.winRect(4) - 70;                    % Posição Y (Canto inferior)
    skipRect(3) = skipRect(1) + 140;                      % Largura do bounding box do botão
    skipRect(4) = skipRect(2) + 40;                       % Altura do bounding box do botão
    
    % Seta invertida no eixo X para apontar para a ESQUERDA (Voltar)
    arrow = [25 -8; -20 -8; -25 0; -20 8; 25 8; 20 0];
    scale = 3.5;
    arrowCoords = arrow * scale + [skipRect(1) + 30, skipRect(2) + 20];

    parentDir = prm.filePath;
    options = mainOptions;
    iconsDir = cellfun(@(x) fullfile(parentDir, prm.imgFolder, [x prm.imgExtension]), options, 'UniformOutput', false);
    clear parentDir

    iconsTex = getMenuTex(dpP.window, iconsDir, drP.black);
    
    leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4};
    
    L = numel(options);

    %% Retângulos dos botões
    gap  = 80;
    btnW = min(300, (dpP.winRect(3)-(L+1)*gap)/L); btnH = 300;
    totalW = L*btnW + (L-1)*gap;
    leftX0 = cX - totalW/2;
    leftX  = leftX0 + (0:L-1) * (btnW + gap);

    btnRects = zeros(4, L);
    btnRects(1,:) = leftX; btnRects(2,:) = cY - btnH/2;
    btnRects(3,:) = leftX + btnW; btnRects(4,:) = cY + btnH/2;
    
    prm.upFrac = 2/3;
    RectSplit = btnRects(2,:)+btnH*prm.upFrac;
    upRects = btnRects; upRects(4,:) = RectSplit; downRects = btnRects; downRects(2,:) = RectSplit;
    uRH = RectHeight(upRects(:,1)'); uRW = RectWidth(upRects(:,1)');
    uRHW = min(uRH, uRW);
    for i=1:L
        [cXaux, cYaux] = RectCenter(upRects(:,i));
        upRects(:,i) = CenterRectOnPoint([0 0, uRHW, uRHW]*prm.iconScaleFactor, cXaux, cYaux);
    end
    clear cXaux cYaux

    %% Retângulo de info adicional
    detailRect = CenterRectOnPoint([0 0 2*prm.btnW (1.2/3)*prm.btnH], cX, 1.75*cY);
    drawBox = false;
    
    %% Estados
    if (L+1)/2 == floor((L+1)/2)
        firstSelectL = (L+1)/2;
        firstSelectR = (L+1)/2;
    else
        firstSelectL = floor((L+1)/2);
        firstSelectR = ceil((L+1)/2);
    end
    selected = 0;
    [prevMx, prevMy, prevBut] = GetMouse(dpP.window);
    prevMouse = [prevMx, prevMy, prevBut];
    
    isMouseMostRecent = false;
    keyWasDown = false;
    lastActionTime = 0; lastRepeatTime = -Inf;
    doAction = false;
    mode = '';
    quitSes = 0;

    while true



        if refreshLayout
            if strcmp(currentScreen, 'main')
                options = mainOptions; optionsName = mainOptionsName; optionsMsg = mainOptionsMsg;
                titleText = 'MENU';
            else
                options = trainOptions; optionsName = trainOptionsName; optionsMsg = trainOptionsMsg;
                titleText = 'TREINO';
            end
            
            % Fecha texturas anteriores para evitar estouro de memória de vídeo
            if exist('iconsTex', 'var'); for i=1:numel(iconsTex); Screen('Close', iconsTex(i)); end; end
            
            L = numel(options);
            btnW = min(300, (dpP.winRect(3)-(L+1)*gap)/L);
            totalW = L*btnW + (L-1)*gap;
            leftX0 = cX - totalW/2;
            leftX  = leftX0 + (0:L-1) * (btnW + gap);
            btnRects = zeros(4, L);
            btnRects(1,:) = leftX; btnRects(2,:) = cY - btnH/2;
            btnRects(3,:) = leftX + btnW; btnRects(4,:) = cY + btnH/2;
            
            prm.upFrac = 2/3;
            RectSplit = btnRects(2,:)+btnH*prm.upFrac;
            upRects = btnRects; upRects(4,:) = RectSplit; downRects = btnRects; downRects(2,:) = RectSplit;
            uRH = RectHeight(upRects(:,1)'); uRW = RectWidth(upRects(:,1)');
            uRHW = min(uRH, uRW);
            for i=1:L
                [cXaux, cYaux] = RectCenter(upRects(:,i));
                upRects(:,i) = CenterRectOnPoint([0 0, uRHW, uRHW]*prm.iconScaleFactor, cXaux, cYaux);
            end
            
            % Carrega novos ícones baseados na pasta de imagens
            iconsDir = cellfun(@(x) fullfile(prm.filePath, prm.imgFolder, [x prm.imgExtension]), options, 'UniformOutput', false);
            iconsTex = getMenuTex(dpP.window, iconsDir, drP.black);
            
            % Reseta seleções
            selected = 0; backSelected = false; refreshLayout = false;
        end

        % Lê posição do mouse
        [mx, my, buttons] = GetMouse(dpP.window);
        mouseMoved = norm([mx my] - prevMouse(1:2)) > 1;
        prevMouse = [mx my buttons];
        cursorVisible = false;
    
        if mouseMoved || isMouseMostRecent
            cursorVisible = true;
            isMouseMostRecent = true;
            selected = 0;
            for i=1:L
                if IsInRect(mx, my, btnRects(:,i))
                    selected = i;
                    break;
                end
            end
            % Verifica botão Voltar (Apenas na tela de treino)
            if strcmp(currentScreen, 'training') && IsInRect(mx, my, skipRect)
                backSelected = true;
                selected = 0;
            end
        end
    
        %% Desenho

        Screen('TextSize', dpP.window, prm.textSizeTitle);
        Screen('DrawText', dpP.window, titleText, titleMargin, titleMargin,  drP.blackGrey);
        Screen('TextSize', dpP.window, prm.textSizeNormal);
        if selected == 0, drawBox = false; end

        % Desenha os retângulos (e seta) e textos e texturas
        Screen('FillRect', dpP.window, drP.blue, btnRects);

        Screen('DrawTextures', dpP.window, iconsTex, [], upRects);

        Screen('TextStyle', dpP.window, 1);
        for i=1:L
            DrawFormattedText(dpP.window, optionsName{i}, 'center', 'center',  drP.black, [], [], [], [], [], downRects(:,i)');
        end
        Screen('TextSize', dpP.window, prm.textSizeNormalish);

        % Se houver algum retângulo selecionado, adiciona contorno
        if selected ~= 0
            drawBox = true;
            msg = optionsMsg{selected};
            Screen('FrameRect', dpP.window,  drP.darkBlue, btnRects(:, selected), prm.pW2);
        end
        Screen('TextStyle', dpP.window, 0);

        if drawBox
            Screen('FillRect', dpP.window, drP.whiteGrey, detailRect);
            Screen('TextSize', dpP.window, prm.textSizeNormalish);
            DrawFormattedText(dpP.window, msg, 'center', 'center',  drP.blackGrey, [], [], [], [], [], detailRect);
            Screen('TextSize', dpP.window, prm.textSizeNormal);
        end

        % Desenha o cursor antes de atualizar a tela
        if cursorVisible
            Screen('FillOval', dpP.window, prm.cursorColor, ...
                [mx-prm.cursorRadius_px my-prm.cursorRadius_px mx+prm.cursorRadius_px my+prm.cursorRadius_px]);
        end

        % --- DESENHA O BOTÃO VOLTAR SE ESTIVER NO SUBMENU ---
        if strcmp(currentScreen, 'training')
            if backSelected
                Screen('FillPoly', dpP.window, drP.darkBlue, arrowCoords);
                Screen('DrawText', dpP.window, skipText, skipRect(1)+60, skipRect(2)+8, drP.darkBlue);
            else
                Screen('FillPoly', dpP.window, drP.blue, arrowCoords);
                Screen('DrawText', dpP.window, skipText, skipRect(1)+60, skipRect(2)+8, drP.black);
            end
        end
    
        Screen('Flip', dpP.window);
    
        % Lê inputs do teclado e clique do mouse
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown || any(buttons)
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
            else
                doAction = true;
            end

            if doAction
                doAction = false;
                if keyCode(leftKey)
                    if selected == 0
                        selected = firstSelectL;
                    else
                        selected = max(1, selected - 1);
                    end
%                     cursorVisible = false;
                    isMouseMostRecent = false;
                elseif keyCode(rightKey)
                    if selected == 0
                        selected = firstSelectR;
                    else
                        selected = min(L, selected + 1);
                    end
%                     cursorVisible = false;
                    isMouseMostRecent = false;
                elseif keyCode(spaceKey) || any(buttons)
                    KbReleaseWait;
                    while any(buttons); [~, ~, buttons] = GetMouse(dpP.window); WaitSecs(0.001); end
                    
                    % Ação 1: Clicou no botão Voltar
                    if strcmp(currentScreen, 'training') && backSelected
                        currentScreen = 'main';
                        refreshLayout = true;
                        continue;
                    end
                    
                    % Ação 2: Selecionou um botão padrão do menu
                    if selected ~= 0
                        mode = options{selected};
                        
                        if strcmp(currentScreen, 'main')
                            % Lógica da Tela Principal
                            if strcmp(mode, 'training')
                                currentScreen = 'training';
                                refreshLayout = true;
                                continue;
                            elseif strcmp(mode, 'staircase')
                                for i=1:L, Screen('Close', iconsTex(i)); end
                                [tkP, taskState] = runStaircase(tkP, dpP, drP, txP, prm, debug, taskState);
                            elseif strcmp(mode, 'experiment')
                                taskState(2,1) = 1;
                                for i=1:L, Screen('Close', iconsTex(i)); end
                                [tkP, taskState, results] = runForaging1(tkP, dpP, drP, txP, prm, debug, mode, taskState);
                            end
                        else
                            % Lógica da Tela de Sub-Treinos
                            taskState(1,1) = 1;
                            for i=1:L, Screen('Close', iconsTex(i)); end
                            [~, taskState, ~] = runForaging1(tkP, dpP, drP, txP, prm, debug, mode, taskState);
                        end
                        
                        % Rotina pós-retorno de qualquer tarefa
                        if taskState(2,2) == 1
                            break; % Encerra se o experimento acabou
                        else
                            refreshLayout = true; % Força recarregamento da tela atual
                            SetMouse(cX, cY, dpP.window);
                        end
                    end
                elseif keyCode(escapeKey)
                    KbReleaseWait;
                    decision = pauseScreen(tkP, dpP, drP, prm, 'menu');
                    if strcmp(decision, 'quit')
                        Screen('Flip', dpP.window);
                        quitSes = 1;
                        break;
                    end
                end
            end
        else
            lastRepeatTime = -Inf;
        end
        keyWasDown = keyIsDown;
    end

    if debug ~= 0, taskState(:) = 0; end
    if ~exist('results', 'var'), results = []; end
    foragingSave(taskState, quitSes, prm, dpP, drP, tkP, txP, results);
end
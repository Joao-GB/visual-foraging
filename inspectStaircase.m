function inspectStaircase(tkP, dpP, drP, prm, RF, thrs)

    screenW = dpP.winRect(3); screenH = dpP.winRect(4);
    targetW = screenW / 2;     targetH = screenH / 2;
    
    % Retângulo de destino
    targetRect = [dpP.winCenter(1) - targetW/2, ...
                  dpP.winCenter(2) - targetH/2, ...
                  dpP.winCenter(1) + targetW/2, ...
                  dpP.winCenter(2) + targetH/2];
    
    B = numel(RF);
    texArray = zeros(1, B+1); % Array para armazenar todas as texturas

    for b = 1:B
        % Cria uma figura invisível
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
        trialNum = 1:length(RF(b).x);
        presentedSigma = -RF(b).x;
        hold on;
        
        % Trajetória
        plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
        
        % Trials corretos
        plot(trialNum(RF(b).response == 1), presentedSigma(RF(b).response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        
        % Trials incorretos
        plot(trialNum(RF(b).response == 0), presentedSigma(RF(b).response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
        
        % Estimativa final do limiar
        yline(thrs(b), '--k', 'LineWidth', 2);
        
        xlabel('Trial'); ylabel('Sigma');
        title(sprintf('Staircase: %s', prm.allOriName{b}))
        grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([1 length(trialNum)]);

        figFrame = getframe(hFig);
        texArray(b) = Screen('MakeTexture', dpP.window, figFrame.cdata);

        baseNameIndiv = sprintf('sc_%s_cond%d', tkP.sesSub, b);
        filenameIndiv = fileVersion(baseNameIndiv, '.png');
        saveas(hFig, filenameIndiv);
        
        close(hFig);
    end

    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
    for b = 1:B
        trialNum = 1:length(RF(b).x);
        presentedSigma = -RF(b).x;
        subplot(1, B, b); hold on;
        
        plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
        plot(trialNum(RF(b).response == 1), presentedSigma(RF(b).response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        plot(trialNum(RF(b).response == 0), presentedSigma(RF(b).response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
        yline(thrs(b), '--k', 'LineWidth', 2);
        
        xlabel('Trial'); title(sprintf('Staircase: %s', prm.allOriName{b}))
        if b == 1, ylabel('Sigma'); end
        grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([1 length(trialNum)]);
    end
    figFrame = getframe(hFig);
    texArray(B+1) = Screen('MakeTexture', dpP.window, figFrame.cdata);
    close(hFig);

    % --- 3. LOOP INTERATIVO DE NAVEGAÇÃO ---
    leftKey   = KbName('LeftArrow'); 
    rightKey  = KbName('RightArrow'); 
    escapeKey = KbName('ESCAPE');
    
    currentView = 1; % Começa na primeira figura individual
    KbReleaseWait;
    
    while true
        % Desenha a textura do estado atual
        Screen('DrawTexture', dpP.window, texArray(currentView), [], targetRect);
        Screen('TextSize', dpP.window, prm.textSizeNormalish);
        
        % Determina o texto explicativo baseado na tela atual
        if currentView <= B
            promptText = 'Pressione uma das setas para avançar pelas figuras';
        else
            promptText = 'Pressione uma das setas para avançar pelas figuras \nou pressione Esc para sair';
        end
        
        textY = targetRect(4) + 50; 
        DrawFormattedText(dpP.window, promptText, 'center', textY, drP.black);
        Screen('Flip', dpP.window);
        
        % Captura comandos do teclado
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(leftKey)
                currentView = max(1, currentView - 1); % Recua sem wrap-around
                KbReleaseWait;
            elseif keyCode(rightKey)
                currentView = min(B+1, currentView + 1); % Avança até a última tela (combinada)
                KbReleaseWait;
            elseif keyCode(escapeKey)
                KbReleaseWait;
                break; % Sai do loop e encerra a inspeção
            end
        end
        WaitSecs(0.01);
    end
    
    % --- 4. LIMPEZA DA MEMÓRIA DE VÍDEO ---
    for i = 1:numel(texArray)
        Screen('Close', texArray(i));
    end
    Screen('Flip', dpP.window);
end
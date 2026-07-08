function inspectStaircase(tkP, dpP, drP, prm, PM, thrs, newThrs, ori)

    screenW = dpP.winRect(3); screenH = dpP.winRect(4);
    targetW = screenW / 2;     targetH = screenH / 2;
    
    % Retângulo de destino
    targetRect = [dpP.winCenter(1) - targetW/2, ...
                  dpP.winCenter(2) - targetH/2, ...
                  dpP.winCenter(1) + targetW/2, ...
                  dpP.winCenter(2) + targetH/2];
    
    B = numel(PM);
    texArray = zeros(1, B + 2 + tkP.stairBurnIn); % Array para armazenar todas as texturas.
                                                  % Incrementa em 1 apenas se houver burn-in

    %% Texturas individuais para cada staircase
    hFigs = gobjects(B,1);
    for b = 1:B
        % Cria uma figura invisível
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
        hFigs(b) = hFig;

        trialNum = 1:length(PM(b).x);
        presentedSigma = -PM(b).x;
        hold on;
        
        % Trajetória
        plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
        
        % Trials corretos
        plot(trialNum(PM(b).response == 1), presentedSigma(PM(b).response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        
        % Trials incorretos
        plot(trialNum(PM(b).response == 0), presentedSigma(PM(b).response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
        
        % Estimativa final do limiar
        if thrs(b) >= newThrs(b)
            yline(thrs(b), '--k', sprintf('75%%: %.2f', thrs(b)), 'LineWidth', 1.5);
            yline(newThrs(b), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(b)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 3);
        else
            yline(thrs(b), '--k', sprintf('75%%: %.2f', thrs(b)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
            yline(newThrs(b), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(b)), 'LineWidth', 3);
        end

        if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev)
            yline(tkP.stairPrev(b).aSigma, '-', 'Prev aSigma', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
        end

        xlabel('Trial'); ylabel('Sigma');
        title(sprintf('Staircase: %s', prm.allOriName{prm.allOriMap(ori(b))}))
        grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([1 length(trialNum)]);

        figFrame = getframe(hFig);
        texArray(b) = Screen('MakeTexture', dpP.window, figFrame.cdata);
    end

    savefig(hFigs, prm.tempFig); close(hFigs); clear hFigs;
    h = openfig(prm.tempFig, 'invisible');
    set(h, 'Visible', 'on');
    savefig(h, prm.tempFig); close(h); clear h;

    
    %% Textura com staircases como subplots
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
    for b = 1:B
        trialNum = 1:length(PM(b).x);
        presentedSigma = -PM(b).x;
        subplot(1, B, b); hold on;
        
        plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
        plot(trialNum(PM(b).response == 1), presentedSigma(PM(b).response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        plot(trialNum(PM(b).response == 0), presentedSigma(PM(b).response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

        if thrs(b) >= newThrs(b)
            yline(thrs(b), '--k', sprintf('75%%: %.2f', thrs(b)), 'LineWidth', 1.5);
            yline(newThrs(b), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(b)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 3);
        else
            yline(thrs(b), '--k', sprintf('75%%: %.2f', thrs(b)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
            yline(newThrs(b), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(b)), 'LineWidth', 3);
        end

        if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev)
            yline(tkP.stairPrev(b).aSigma, '-', 'Prev aSigma', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
        end
        
        xlabel('Trial'); title(sprintf('Staircase: %s', prm.allOriName{prm.allOriMap(ori(b))}))
        if b == 1, ylabel('Sigma'); end
        grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([1 length(trialNum)]);
    end
    figFrame = getframe(hFig);
    texArray(B+1) = Screen('MakeTexture', dpP.window, figFrame.cdata);
    close(hFig);

    %% Textura com staircases como subplots, sem burn-in
    if tkP.stairBurnIn
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
        for b = 1:B
            usefulTrials = ((prm.nStimsStair-1)*prm.burninTrials+1):length(PM(b).x);
            trialNum = usefulTrials;
            presentedSigma = -PM(b).x(usefulTrials);
            response = PM(b).response(usefulTrials);
            subplot(1, B, b); hold on;
            
            plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
            plot(trialNum(response == 1), presentedSigma(response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
            plot(trialNum(response == 0), presentedSigma(response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
            if thrs(b) >= newThrs(b)
                yline(thrs(b), '--k', sprintf('75%%: %.2f', thrs(b)), 'LineWidth', 1.5);
                yline(newThrs(b), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(b)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 3);
            else
                yline(thrs(b), '--k', sprintf('75%%: %.2f', thrs(b)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
                yline(newThrs(b), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(b)), 'LineWidth', 3);
            end
            xlabel('Trial'); title(sprintf('Staircase: %s', prm.allOriName{prm.allOriMap(ori(b))}))
            if b == 1, ylabel('Sigma'); end
            grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([trialNum(1) trialNum(end)]);
        end
        figFrame = getframe(hFig);
        texArray(B+2) = Screen('MakeTexture', dpP.window, figFrame.cdata);
        close(hFig);
    end

    %% Textura com Curvas Psicométricas como subplots
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
    for b = 1:B
        subplot(1, B, b); hold on;
        
        % Extrai dados do staircase
        presentedSigma = -PM(b).x; 
        
        % Agrupa trials por intensidade para plotar os círculos proporcionais
        [SL, NP, OON] = PAL_PFML_GroupTrialsbyX(presentedSigma, PM(b).response, ones(size(PM(b).response)));
        for SR = 1:length(SL(OON~=0))
            plot(SL(SR), NP(SR)/OON(SR), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 20*sqrt(OON(SR)/sum(OON)));
        end
        
        % Plota a curva ajustada
        stimRange = linspace(prm.sigmaMin, prm.sigmaMax, 100);
        alphaEst = PM(b).threshold(end);
        betaEst  = PM(b).slope(end);
        plot(stimRange, PAL_CumulativeNormal([alphaEst, betaEst, PM(b).gamma(end), PM(b).lapse(end)], stimRange), 'k-', 'LineWidth', 2);
        
        % Insere o texto com os parâmetros
        txt = sprintf('Curr:\nThrs: %.2f\nSlope: %.2f', alphaEst, betaEst);
        if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev)
            txt = sprintf('%s\n\nPrev:\nThrs: %.2f\nSlope: %.2f', txt, ...
                          tkP.stairPrev(b).threshold(end), tkP.stairPrev(b).slope(end));
        end
        % Posiciona o texto no canto superior esquerdo do plot
        text(prm.sigmaMin + (prm.sigmaMax-prm.sigmaMin)*0.05, 0.95, txt, 'FontSize', 9, 'VerticalAlignment', 'top');
        
        xlabel('Sigma'); title(sprintf('Curve: %s', prm.allOriName{prm.allOriMap(ori(b))}));
        if b == 1, ylabel('Proportion Correct'); end
        grid on; ylim([0 1.05]); xlim([prm.sigmaMin prm.sigmaMax]);
    end
    figFrame = getframe(hFig);
    texArray(B + 2 + tkP.stairBurnIn) = Screen('MakeTexture', dpP.window, figFrame.cdata);
    close(hFig);

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
                currentView = min(length(texArray), currentView + 1); % Avança até a última tela (combinada)
                KbReleaseWait;
            elseif keyCode(escapeKey)
                KbReleaseWait;
                break; % Sai do loop e encerra a inspeção
            end
        end
        WaitSecs(0.01);
    end
    
    % Limpeza de memória: exclui texturas
    for i = 1:numel(texArray)
        Screen('Close', texArray(i));
    end
    Screen('Flip', dpP.window);
end
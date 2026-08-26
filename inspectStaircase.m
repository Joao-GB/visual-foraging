function newNewThrs = inspectStaircase(tkP, dpP, drP, prm, PM, ~, newThrs, ori)
% Os thrs, newThrs e tkP.stairPrev estão na ordem canônica de estímulos, 
% mas o PM não, por isso a estratégia do k
    oldNewThrs = aSigmaFromStair(newThrs, prm, 1);
    screenW = dpP.winRect(3); screenH = dpP.winRect(4);
    targetW = screenW / 2;    targetH = screenH / 2;
    
    targetRect = [dpP.winCenter(1) - targetW/2, ...
                  dpP.winCenter(2) - targetH/2, ...
                  dpP.winCenter(1) + targetW/2, ...
                  dpP.winCenter(2) + targetH/2];
    
    B = numel(PM); % Número de blocos
    
    % Total de texturas que serão geradas (B blocos + 1 curva + 1 opcional burn-in)
    nTextures = B + tkP.stairBurnIn; % B + 1 + tkP.stairBurnIn;
    texArray = zeros(1, nTextures); 
    
    % Apenas B+1 figuras serão salvas no .fig (B blocos + 1 curva conjunta)
    hFigs = gobjects(B,1);
    
    %% Texturas 1 a B: Staircases individuais com histórico completo
    for b = 1:B
        tgtOri = PM(b).tgtOri;
        if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev), k = find([tkP.stairPrev.tgtOri] == tgtOri,1); end
        l = prm.allOriMap(tgtOri);
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
        hFigs(b) = hFig;
        
        trialNum = 1:length(PM(b).x);
        presentedSigma = -PM(b).x;
        hold on;
        
        % Trajetória
        plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
        
        % Trials corretos e incorretos
        plot(trialNum(PM(b).response == 1), presentedSigma(PM(b).response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        plot(trialNum(PM(b).response == 0), presentedSigma(PM(b).response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
        
        % Estimativa final do limiar
        yline(newThrs(l), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(l)), 'LineWidth', 3);
        if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev)
            yline(tkP.stairPrev(k).aSigma, '-', 'Prev aSigma', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
        end
        xlabel('Trial'); ylabel('Sigma');
        title(sprintf('Staircase: %s', prm.allOriName{prm.allOriMap(ori(b))}))
        grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([1 length(trialNum)]);
        
        figFrame = getframe(hFig);
        texArray(b) = Screen('MakeTexture', dpP.window, figFrame.cdata);
    end
    
    % Variável para garantir que as texturas seguintes fiquem nas posições certas
    texIdx = B + 1;
    
    %% Textura: Staircase sem burn-in (se houver)
    if tkP.stairBurnIn
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
        
        for b = 1:B
            tgtOri = PM(b).tgtOri;
            if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev), k = find([tkP.stairPrev.tgtOri] == tgtOri,1); end
            usefulTrials = ((prm.nStimsStair-1)*prm.burninTrials+1):length(PM(b).x);
            trialNum = usefulTrials;
            presentedSigma = -PM(b).x(usefulTrials);
            response = PM(b).response(usefulTrials);
            subplot(1, B, b); hold on;
            
            plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
            plot(trialNum(response == 1), presentedSigma(response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
            plot(trialNum(response == 0), presentedSigma(response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
    
            yline(newThrs(k), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(k)), 'LineWidth', 3);
            xlabel('Trial'); title(sprintf('Staircase: %s', prm.allOriName{prm.allOriMap(ori(b))}))
            if b == 1, ylabel('Sigma'); end
            grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([trialNum(1) trialNum(end)]);
        end
        figFrame = getframe(hFig);
        texArray(texIdx) = Screen('MakeTexture', dpP.window, figFrame.cdata);
%         texIdx = texIdx + 1;
        close(hFig); % Fecha sem salvar no array hFigs
    end
    
    %% Save & Display Workflow
    savefig(hFigs, prm.tempFig); close(hFigs); clear hFigs;
    h = openfig(prm.tempFig, 'invisible');
    set(h, 'Visible', 'on');
    savefig(h, prm.tempFig); close(h); clear h;
    
    %% Interação via PTB
    leftKey   = KbName('LeftArrow');
    rightKey  = KbName('RightArrow');
    upKey   = KbName('UpArrow');
    downKey  = KbName('DownArrow');
    escapeKey = KbName('ESCAPE');
    
    currentView = 1; 
    KbReleaseWait;

    newNewThrs = oldNewThrs;
    textColor1   = repmat(drP.black, 1, 3);
    while true
        % Draw the current figure
        Screen('DrawTexture', dpP.window, texArray(currentView), [], targetRect);
    
        if currentView >= length(texArray)
    
            % Small box dimensions
            boxWidth  = 80;
            boxHeight = 35;
            boxMargin = 20;  % distance between figure and threshold box
    
            % Position: upper-right of the figure, just outside it
            boxLeft = targetRect(3) + boxMargin;
            boxTop  = targetRect(2);
    
            boxRect = [boxLeft, ...
                       boxTop, ...
                       boxLeft + boxWidth, ...
                       boxTop + boxHeight];
    
            % White box
            Screen('FillRect', dpP.window, [255 255 255], boxRect);
    
            % Threshold text inside the box
            Screen('TextSize', dpP.window, prm.textSizeNormalish);
    
            thresholdText = num2str(newNewThrs);
            DrawFormattedText(dpP.window, thresholdText, ...
                'center', ...
                'center', ...
                drP.black, [], [], [], [], [], boxRect);

    
            % -----------------------------------------------------
            % Up triangle above the box
            % -----------------------------------------------------
            triangleWidth  = 14;
            triangleHeight = 10;
    
            triangleCenterX = boxLeft + boxWidth/2;
    
            upTriangle = [ ...
                triangleCenterX,                    boxTop - triangleHeight; ...
                triangleCenterX - triangleWidth/2, boxTop; ...
                triangleCenterX + triangleWidth/2, boxTop ...
            ]';
    
            Screen('FillPoly', dpP.window, textColor1, upTriangle');
    
            % -----------------------------------------------------
            % Down triangle below the box
            % -----------------------------------------------------
            downTriangle = [ ...
                triangleCenterX,                    boxTop + boxHeight + triangleHeight; ...
                triangleCenterX - triangleWidth/2, boxTop + boxHeight; ...
                triangleCenterX + triangleWidth/2, boxTop + boxHeight ...
            ]';
    
            Screen('FillPoly', dpP.window, textColor1, downTriangle');
    
            % -----------------------------------------------------
            % "(recommended)" if threshold equals oldNewThrs
            % -----------------------------------------------------
            if newNewThrs == oldNewThrs
    
                Screen('TextSize', dpP.window, prm.textSizeSmall);
    
                recommendedText = '(recomendado)';
    
                % Put it immediately to the right of the box
                recommendedX = boxRect(3) + 8;
                recommendedY = boxTop + boxHeight/2 - prm.textSizeSmall/2 + 6;
                
                DrawFormattedText(dpP.window, recommendedText, ...
                    recommendedX, recommendedY, drP.black);

            end
        end
    
        % ---------------------------------------------------------
        % Prompt
        % ---------------------------------------------------------
        Screen('TextSize', dpP.window, prm.textSizeNormalish);
    
        if currentView < length(texArray)
            promptText = 'Pressione uma das setas para avançar pelas figuras';
        else
            promptText = ['Pressione uma das setas para avançar pelas figuras ' ...
                          '\nou pressione Esc para sair'];
        end
    
        textY = targetRect(4) + 50;
        DrawFormattedText(dpP.window, promptText, ...
            'center', textY, drP.black);
    
        Screen('Flip', dpP.window);
    
        % ---------------------------------------------------------
        % Keyboard
        % ---------------------------------------------------------
        [keyIsDown, ~, keyCode] = KbCheck;
    
        if keyIsDown
    
            if keyCode(leftKey)
    
                currentView = max(1, currentView - 1);
                KbReleaseWait;
    
            elseif keyCode(rightKey)
    
                currentView = min(length(texArray), currentView + 1);
                KbReleaseWait;
    
            elseif keyCode(escapeKey)
    
                KbReleaseWait;
                break;
    
            % -----------------------------------------------------
            % Change threshold ONLY on the last figure
            % -----------------------------------------------------
            elseif currentView >= length(texArray) && keyCode(upKey)
    
                newNewThrs = min(prm.sigmaMax, ...
                                 newNewThrs + prm.sigmaRem);
                KbReleaseWait;
    
            elseif currentView >= length(texArray) && keyCode(downKey)
    
                newNewThrs = max(prm.sigmaMin, ...
                                 newNewThrs - prm.sigmaRem);
                KbReleaseWait;
            end
        end
    
        WaitSecs(0.01);
    end
    
    % Exclui texturas
    for i = 1:numel(texArray)
        Screen('Close', texArray(i));
    end
    Screen('Flip', dpP.window);
end
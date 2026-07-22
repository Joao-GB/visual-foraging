function inspectStaircase(tkP, dpP, drP, prm, PM, thrs, newThrs, ori)
% Os thrs, newThrs e tkP.stairPrev estão na ordem canônica de estímulos, 
% mas o PM não, por isso a estratégia do k
    screenW = dpP.winRect(3); screenH = dpP.winRect(4);
    targetW = screenW / 2;    targetH = screenH / 2;
    
    targetRect = [dpP.winCenter(1) - targetW/2, ...
                  dpP.winCenter(2) - targetH/2, ...
                  dpP.winCenter(1) + targetW/2, ...
                  dpP.winCenter(2) + targetH/2];
    
    B = numel(PM); % Número de blocos
    
    % Total de texturas que serão geradas (B blocos + 1 curva + 1 opcional burn-in)
    nTextures = B + 1 + tkP.stairBurnIn;
    texArray = zeros(1, nTextures); 
    
    % Apenas B+1 figuras serão salvas no .fig (B blocos + 1 curva conjunta)
    hFigs = gobjects(B+1,1);
    
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
        if thrs(l) >= newThrs(l)
            yline(thrs(l), '--k', sprintf('75%%: %.2f', thrs(l)), 'LineWidth', 1.5);
            yline(newThrs(l), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(l)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 3);
        else
            yline(thrs(l), '--k', sprintf('75%%: %.2f', thrs(l)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
            yline(newThrs(l), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(l)), 'LineWidth', 3);
        end
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
    
            if thrs(k) >= newThrs(k)
                yline(thrs(k), '--k', sprintf('75%%: %.2f', thrs(k)), 'LineWidth', 1.5);
                yline(newThrs(k), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(k)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 3);
            else
                yline(thrs(k), '--k', sprintf('75%%: %.2f', thrs(k)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
                yline(newThrs(k), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(k)), 'LineWidth', 3);
            end
            xlabel('Trial'); title(sprintf('Staircase: %s', prm.allOriName{prm.allOriMap(ori(b))}))
            if b == 1, ylabel('Sigma'); end
            grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([trialNum(1) trialNum(end)]);
        end
        figFrame = getframe(hFig);
        texArray(texIdx) = Screen('MakeTexture', dpP.window, figFrame.cdata);
        texIdx = texIdx + 1;
        close(hFig); % Fecha sem salvar no array hFigs
    end
    
    %% Textura: Curvas psicométricas para o histórico completo do staircase
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
    hFigs(B+1) = hFig; % Guarda na última posição válida para salvar
    
    for b = 1:B
        tgtOri = PM(b).tgtOri;
        if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev), k = find([tkP.stairPrev.tgtOri] == tgtOri,1); end
        subplot(1, B, b); hold on;
    
        nBins = 15;
        binEdges = linspace(prm.sigmaMin, prm.sigmaMax, nBins + 1);
        binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
        rawSigma = -PM(b).x(1:end-1); 
        binIdx = discretize(rawSigma, binEdges);
        binIdx(isnan(binIdx)) = 1; 
        binnedSigma = binCenters(binIdx);
        
        [SL, NP, OON] = PAL_PFML_GroupTrialsbyX(binnedSigma, PM(b).response, ones(size(PM(b).response)));
        for SR = 1:length(SL(OON~=0))
            plot(SL(SR), NP(SR)/OON(SR), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 20*sqrt(OON(SR)/sum(OON)));
        end
        
        stimRange = linspace(-prm.sigmaMax, -prm.sigmaMin, 100);
        alphaEst = PM(b).threshold(end);
        betaEst  = PM(b).slope(end);
        plot(-stimRange, PAL_CumulativeNormal([alphaEst, betaEst, PM(b).guess(end), PM(b).lapse(end)], stimRange), 'k-', 'LineWidth', 2);
        set(gca, 'XDir', 'reverse');
        
        txt = sprintf('Atual:\nLimiar: %.2f\nIncl.: %.2f', -alphaEst, betaEst);
        if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev)
            txt = sprintf('%s\n\nPrev.:\nLimiar: %.2f\nIncl.: %.2f', txt, ...
                          -tkP.stairPrev(k).threshold(end), tkP.stairPrev(k).slope(end));
        end
        text(prm.sigmaMin + (prm.sigmaMax-prm.sigmaMin)*0.05, 0.95, txt, 'FontSize', 9, 'VerticalAlignment', 'top');
        
        xlabel('Sigma'); title(sprintf('Curve: %s', prm.allOriName{prm.allOriMap(ori(b))}));
        if b == 1, ylabel('Proporção corretos'); end
        grid on; ylim([0 1.05]); xlim([prm.sigmaMin prm.sigmaMax]);
    end
    figFrame = getframe(hFig);
    texArray(texIdx) = Screen('MakeTexture', dpP.window, figFrame.cdata);
    
    %% Save & Display Workflow
    savefig(hFigs, prm.tempFig); close(hFigs); clear hFigs;
    h = openfig(prm.tempFig, 'invisible');
    set(h, 'Visible', 'on');
    savefig(h, prm.tempFig); close(h); clear h;
    
    %% Interação via PTB
    leftKey   = KbName('LeftArrow'); 
    rightKey  = KbName('RightArrow'); 
    escapeKey = KbName('ESCAPE');
    
    currentView = 1; 
    KbReleaseWait;
    
    while true
        Screen('DrawTexture', dpP.window, texArray(currentView), [], targetRect);
        Screen('TextSize', dpP.window, prm.textSizeNormalish);
        
        if currentView < length(texArray)
            promptText = 'Pressione uma das setas para avançar pelas figuras';
        else
            promptText = 'Pressione uma das setas para avançar pelas figuras \nou pressione Esc para sair';
        end
        
        textY = targetRect(4) + 50; 
        DrawFormattedText(dpP.window, promptText, 'center', textY, drP.black);
        Screen('Flip', dpP.window);
        
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
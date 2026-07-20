function inspectStaircase(tkP, dpP, drP, prm, PM, thrs, newThrs, ori)
    screenW = dpP.winRect(3); screenH = dpP.winRect(4);
    targetW = screenW / 2;     targetH = screenH / 2;
    
    targetRect = [dpP.winCenter(1) - targetW/2, ...
                  dpP.winCenter(2) - targetH/2, ...
                  dpP.winCenter(1) + targetW/2, ...
                  dpP.winCenter(2) + targetH/2];
    
    B = numel(PM); % Número de blocos
    
    % Encontra os trials onde houve mudança de bloco (tiro 1 pois será
    % produzida uma sugestão que não usarei)
    blockTransitions = zeros(1, B-1);
    for b = 1:(B-1)
        blockTransitions(b) = length(PM(b).x) - 1;
    end
    
    % Array para armazenar as texturas: 1 para staircase inteira, 1 para curva psicométrica, +1 se houver burn-in
    nTextures = 2 + tkP.stairBurnIn; 
    texArray = zeros(1, nTextures); 
    hFigs = gobjects(nTextures, 1);
    
    %% Textura 1: Staircase com histórico completo
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
    hFigs(1) = hFig;
    
    trialNum = 1:length(PM(B).x);
    presentedSigma = -PM(B).x;
    hold on;
    
    % Trajetória
    plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
    
    % Trials corretos e incorretos
    plot(trialNum(PM(B).response == 1), presentedSigma(PM(B).response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    plot(trialNum(PM(B).response == 0), presentedSigma(PM(B).response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
    
    % Linhas verticais para indicar as trocas de bloco
    if ~isempty(blockTransitions)
        xline(blockTransitions + 0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end
    
    % Estimativa final do limiar
    if thrs(B) >= newThrs(B)
        yline(thrs(B), '--k', sprintf('75%%: %.2f', thrs(B)), 'LineWidth', 1.5);
        yline(newThrs(B), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(B)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 3);
    else
        yline(thrs(B), '--k', sprintf('75%%: %.2f', thrs(B)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
        yline(newThrs(B), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(B)), 'LineWidth', 3);
    end
    if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev)
        yline(tkP.stairPrev(B).aSigma, '-', 'Prev aSigma', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    end
    xlabel('Trial'); ylabel('Sigma');
    auxS = prm.allOriName{prm.allOriMap(ori(1))};
    for i = 2:B
        auxS = [auxS ', ' prm.allOriName{prm.allOriMap(ori(i))}]; %#ok<AGROW>
    end
    title(sprintf('Staircase: %s (todos blocos)', auxS))
    grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([1 length(trialNum)]);
    
    figFrame = getframe(hFig);
    texArray(1) = Screen('MakeTexture', dpP.window, figFrame.cdata);
    
    %% Textura 2: Staircase sem burn-in (se houver)
    texIdx = 2;
    if tkP.stairBurnIn
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
        hFigs(texIdx) = hFig;
        
        usefulTrials = ((prm.nStimsStair-1)*prm.burninTrials+1):length(PM(B).x);
        trialNum = usefulTrials;
        presentedSigma = -PM(B).x(usefulTrials);
        response = PM(B).response(usefulTrials);
        hold on;
        
        plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
        plot(trialNum(response == 1), presentedSigma(response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        plot(trialNum(response == 0), presentedSigma(response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

        validTransitions = blockTransitions(blockTransitions >= usefulTrials(1));
        if ~isempty(validTransitions)
            xline(validTransitions + 0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        end
        
        if thrs(B) >= newThrs(B)
            yline(thrs(B), '--k', sprintf('75%%: %.2f', thrs(B)), 'LineWidth', 1.5);
            yline(newThrs(B), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(B)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 3);
        else
            yline(thrs(B), '--k', sprintf('75%%: %.2f', thrs(B)), 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
            yline(newThrs(B), '-k', sprintf('%d%%: %.2f', round(prm.stairLevel*100), newThrs(B)), 'LineWidth', 3);
        end
        xlabel('Trial'); ylabel('Sigma');
        title(sprintf('Staircase: %s (sem burn-in)', auxS))
        grid on; ylim([prm.sigmaMin prm.sigmaMax]); xlim([trialNum(1) trialNum(end)]);
        
        figFrame = getframe(hFig);
        texArray(texIdx) = Screen('MakeTexture', dpP.window, figFrame.cdata);
        texIdx = texIdx + 1;
    end
    
    %% Textura 3: Curvas psicométricas para o histórico completo do staircase
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
    hFigs(texIdx) = hFig;
    hold on;
    
    nBins = 15;
    binEdges = linspace(prm.sigmaMin, prm.sigmaMax, nBins + 1);
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
    rawSigma = -PM(B).x(1:end-1); 
    binIdx = discretize(rawSigma, binEdges);
    binIdx(isnan(binIdx)) = 1; 
    binnedSigma = binCenters(binIdx);
    
    [SL, NP, OON] = PAL_PFML_GroupTrialsbyX(binnedSigma, PM(B).response, ones(size(PM(B).response)));
    for SR = 1:length(SL(OON~=0))
        plot(SL(SR), NP(SR)/OON(SR), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 20*sqrt(OON(SR)/sum(OON)));
    end
    
    stimRange = linspace(-prm.sigmaMax, -prm.sigmaMin, 100);
    alphaEst = PM(B).threshold(end);
    betaEst  = PM(B).slope(end);
    plot(-stimRange, PAL_CumulativeNormal([alphaEst, betaEst, PM(B).guess(end), PM(B).lapse(end)], stimRange), 'k-', 'LineWidth', 2);
    set(gca, 'XDir', 'reverse');
    
    txt = sprintf('Atual:\nLimiar: %.2f\nIncl.: %.2f', -alphaEst, betaEst);
    if isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev)
        txt = sprintf('%s\n\nPrev.:\nLimiar: %.2f\nIncl.: %.2f', txt, ...
                      -tkP.stairPrev(B).threshold(end), tkP.stairPrev(B).slope(end));
    end
    text(prm.sigmaMin + (prm.sigmaMax-prm.sigmaMin)*0.05, 0.95, txt, 'FontSize', 9, 'VerticalAlignment', 'top');
    
    xlabel('Sigma'); ylabel('Proporção corretos');
    title(sprintf('Curva: %s', auxS));
    grid on; ylim([0 1.05]); xlim([prm.sigmaMin prm.sigmaMax]);

    figFrame = getframe(hFig);
    texArray(texIdx) = Screen('MakeTexture', dpP.window, figFrame.cdata);
    
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
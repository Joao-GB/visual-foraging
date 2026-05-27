function inspectStaircase(dpP, drP, prm, RF, thrs)

    screenW = dpP.winRect(3); screenH = dpP.winRect(4);
    targetW = screenW / 2;     targetH = screenH / 2;
    
    % Retângulo de destino
    targetRect = [dpP.winCenter(1) - targetW/2, ...
                  dpP.winCenter(2) - targetH/2, ...
                  dpP.winCenter(1) + targetW/2, ...
                  dpP.winCenter(2) + targetH/2];
% Cria uma figura invisível
    % hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 targetW targetH]);
    hFig = figure('Units', 'pixels', 'Position', [0 0 targetW targetH]);
    
    B = numel(RF);
    for b = 1:B
        trialNum = 1:length(RF(b).x);
        presentedSigma = -RF(b).x;
        
        subplot(1,B,b); hold on;
        
        % Trajetória
        plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
        
        % Trials corretos
        plot(trialNum(RF(b).response == 1), presentedSigma(RF(b).response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        
        % Trials incorretos
        plot(trialNum(RF(b).response == 0), presentedSigma(RF(b).response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
        
        % Estimativa final do limiar
        yline(thrs(b), '--k', 'LineWidth', 2);
        
        xlabel('Trial');
        title(sprintf('Staircase: %s', prm.allOriName{b}))
        if b == 1, ylabel('Sigma'); end
        grid on;
        
        ylim([prm.sigmaMin prm.sigmaMax]);
        xlim([1 length(trialNum)]);
    end

    figFrame = getframe(hFig);
    graphMatrix = figFrame.cdata;
    close(hFig);
    graphTex = Screen('MakeTexture', dpP.window, graphMatrix);
    Screen('DrawTexture', dpP.window, graphTex, [], targetRect);

    Screen('TextSize', dpP.window, prm.textSizeNormalish);
    promptText = 'Pressione qualquer tecla para retornar ao menu';

    textY = targetRect(4) + 50; 
    DrawFormattedText(dpP.window, promptText, 'center', textY, drP.black);

    Screen('Flip', dpP.window);

    KbReleaseWait;
    KbWait;

    Screen('Close', graphTex);
    Screen('Flip', dpP.window);
end
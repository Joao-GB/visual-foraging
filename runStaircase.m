function [resultsStair] = runStaircase(tkP, dpP, drP, txP, prm)
%% Pré
% A tela de estímulos deve ser similar à do experimento, quanto à
% quantidade de estímulos e ao início aleatório. 
        Screen('Flip', dpP.window);
        nBlocks = prm.nBlocksStair; nTrials = prm.nTrialsStair;

        if isfield(tkP,'pinkNoiseDur'), prm.pinkNoiseDur = tkP.pinkNoiseDur; end
        
        leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4}; rKey = tkP.keys{5};
        if strcmp(tkP.targetKey, 'right'), targetKey = rightKey; nonTargetKey = leftKey;
        else,                              targetKey = leftKey; nonTargetKey = rightKey;
        end
        % Não me interessam os outputs de distribuições temporais, pois no
        % staircase o tempo não será por quantidade de fixações
        nTrialsBuffered = nTrials + prm.nBufferTrials;
        [nTs, nStims, targetOri, ~, nStimsToReport, ~] = getForagingDistributions1(tkP.nStims, tkP.nMinFix, tkP.nMaxFix, nTrialsBuffered, nBlocks, prm);

        nStimsToReport(:) = 3;
        drP.allColors = drP.white*ones(3, nStims);
        drP.allPW     = prm.pW1*ones(1,nStims);
        
        orientation = zeros(nStims, nTrialsBuffered, nBlocks);
        for b=1:nBlocks
            lowerBound = targetOri(b) + prm.nbhdRadius;
            upperBound = 180 + (targetOri(b) - prm.nbhdRadius);
            orientation(:,:,b) = mod(rand(nStims, nTrialsBuffered)*(upperBound-lowerBound)+lowerBound, 180);
        end

        crossSize_px = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.crossSize_dva);
        fixCenters  = zeros(2, nTrialsBuffered, nBlocks);

        noiseMatrix = zeros(txP.gabor.size_px, (nStims)*txP.gabor.size_px);
        oriPinkMatrix = zeros(txP.gabor.size_px, (nStims)*txP.gabor.size_px);
    
        noiseCenters = [0:txP.gabor.size_px:(nStims-1)*txP.gabor.size_px; zeros(1, nStims)] + txP.gabor.size_px/2;
        baseRect = [0 0 txP.gabor.size_px txP.gabor.size_px];
    
        srcRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', noiseCenters(1,:), noiseCenters(2,:));
    
        minDist_px   = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.minDist_dva+prm.gaborSize_dva);
        minFixDist1 = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.fixROIradius1_dva);
    
        stimCenters = zeros(2, nStims, nTrialsBuffered, nBlocks);
        fixIdx      = zeros(nBlocks, nTrialsBuffered); 
        for b=1:nBlocks
            for i=1:nTrialsBuffered
                [currFixCenter, currStimCenter, ~, fI] = getStimLocations2_1(dpP.winRect(3:4), [dpP.winCenter 1], nStims, minDist_px, txP.gabor.size_px);
                fixCenters(:, i, b)     = currFixCenter;
                stimCenters(:, :, i, b) = currStimCenter;
                fixIdx(b, i) = fI;
                orientation(randperm(nStims, nTs(b, i)), i, b) = targetOri(b);
            end
        end
    
        minFixDist3 = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.fixROIradius3_dva);

%% Cria o objeto que registra todo o staircase
        alphaRange = (-prm.sigmaMax:.1:-prm.sigmaMin);
        PF = @PAL_CumulativeNormal;
        prior = PAL_pdfNormal(alphaRange, -prm.priorMeanStair, prm.priorStdStair);
        RF = PAL_AMRF_setupRF( ...
            'priorAlphaRange', alphaRange, ...
            'prior', prior, ...
            'PF', PF, ...
            'beta', 1, 'gamma', 0.5, 'lambda', 0.02, ...
            'xMin', -prm.sigmaMax, ...
            'xMax', -prm.sigmaMin, ...
            'meanmode', 'mean', ...
            'stopCriterion', 'trials', ...
            'stopRule', 3*nTrials);
        RF(1:nBlocks) = RF(1);
        aSigma = ones(1, nBlocks)*prm.aSigma;

        oriFilter = repmat(txP.oriFilter, [1, 1, nBlocks]);
        OFsize = repmat(txP.OFsize, [1, 1, nBlocks]);
        jitterTimes = rand(nBlocks, nTrialsBuffered)*(prm.maxJitterStair-prm.minJitterStair)+prm.minJitterStair;

%% Tarefa
        fprintf('----Início do staircase----\n')
        try
%% Tela de blocos
            Screen('TextFont', dpP.window, prm.textFont);
            keepGoingBlocks = true;
            wasRecording = false;
            restartBlock = false;
            b = 1;

            trialOrder = zeros(2, nTrials, nBlocks);
            trialFeedback = cell(nBlocks, nTrials);

            Eyelink('Message',prm.msg.on.stc);
            EyelinkDoTrackerSetup(tkP.el);

            while b <= nBlocks && keepGoingBlocks
                trialOrder(:,:,b) = 0;

                overlayRect = CenterRectOnPointd([dpP.winRect(1:3) dpP.winRect(4)/2.5], dpP.winRect(3)/2, dpP.winRect(4)/2);
                overlayColor = repmat(drP.whiteGrey, 1, 3);
                textColor1   = repmat(drP.black, 1, 3);
                textColor2   = repmat(.5*(drP.grey+drP.black), 1, 3);

                repeatMessage = true;
                
                while repeatMessage && keepGoingBlocks
                    Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
                    Screen('DrawTexture', dpP.window, txP.exampleGabor.tex, [], [], targetOri(b), [], [], [], [], []);
                    Screen('DrawTexture', dpP.window, txP.exampleNoise.tex, [], [], [], [], [], [], [], []);
                    Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                    Screen('DrawTextures', dpP.window, txP.exampleBlob.tex, [], [], targetOri(b), [], [], [0 0 0 1]', [], [], txP.exampleBlob.props); 
        
                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    Screen('TextSize', dpP.window, prm.textSizeBigger);
                    countBlockText  = sprintf('Bloco %d de %d', b, nBlocks);
                    DrawFormattedText(dpP.window, countBlockText, 'center', dpP.screenRes.height*.1, drP.black);
        
                    Screen('TextSize', dpP.window, prm.textSizeHuge);
                    targetText_1    = 'Orientação dos alvos:';
                    DrawFormattedText(dpP.window, targetText_1, 'center', dpP.screenRes.height*.2, drP.black);
        
                    targetText_2    = prm.allOriName{prm.allOriMap(targetOri(b))};
                    Screen('TextStyle', dpP.window, 1); Screen('TextSize', dpP.window, prm.textSizeHuger);
                    DrawFormattedText(dpP.window, targetText_2, 'center', dpP.screenRes.height*.3, drP.black);
                                                                    
                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    Screen('TextStyle', dpP.window, 0);
                    Screen('TextSize', dpP.window, prm.textSizeHuge);
                    proceedText = 'Pressione ESPAÇO para prosseguir';
                    DrawFormattedText(dpP.window, proceedText, 'center', dpP.screenRes.height*.75, drP.black);
        
                    Screen('TextSize', dpP.window, prm.textSizeNormal);
        
                    Screen('Flip', dpP.window);
        
                    KbReleaseWait;
                    KbWait;
                    
                    [~,~,keyCode] = KbCheck;
                    
                    if keyCode(spaceKey)
                        KbReleaseWait;
                        repeatMessage = false;
                    
                    elseif keyCode(escapeKey)
                        KbReleaseWait;
                        [keepGoingBlocks, ~, ~, ~] = pauseHandle(keepGoingBlocks, [], [], [], wasRecording, tkP, txP, dpP, drP, prm, 'block', 0, 2, []);
                    end
                end

                Screen('Flip', dpP.window);
                
                if keepGoingBlocks
                    Eyelink('Message', sprintf(prm.msg.on.blk{1}, b, nBlocks));
                end
                i = 1;

                trialQueue = randperm(nTrialsBuffered);
                retryCount = zeros(1, nTrialsBuffered);

                keepGoingTrials = keepGoingBlocks;
%% Loop principal
                while i <= nTrials && keepGoingTrials
                    trialIdxUp = false;
                    idx = trialQueue(i);
                    trialOrder(1, i, b) = idx;
                    restartTrial = false;
                    fprintf('\nIdx Bloco: %d\n# Trial: %d\n Idx Trial: %d\n', b, i, idx)
                    dstRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', stimCenters(1,:, idx, b), stimCenters(2,:, idx, b));
                
                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    Eyelink('Command', 'clear_screen 0');
                    for j=1:size(dstRects, 2)
                        hostRect = dstRects(:,j);
                        cX = mean(hostRect([1 3]));
                        cY = mean(hostRect([2 4]));
                        hostRect = round(hostRect);
                        radius = (hostRect(3)-hostRect(1))/2;
                        theta = orientation(j, idx, b); lineLen = radius * 0.8;
                        dx = lineLen * cosd(theta); dy = lineLen * sind(theta);
                        Eyelink('Command', 'draw_line %d %d %d %d 12', round(cX-dx), round(cY-dy), round(cX+dx), round(cY+dy));
                    end
                    
%% Fase 1: tela de fixação
                    auxNoiseMatrix = butterFilter(pinkNoise(txP.gabor.size_px, (nStims)*txP.gabor.size_px), txP.noiseLoCutFreq, txP.noiseHiCutFreq);
                    for j=1:(nStims)
                        colRange = ((j-1)*txP.gabor.size_px+1):(j*txP.gabor.size_px);
                        aux = auxNoiseMatrix(:,colRange);
                        noiseMatrix(:, colRange) = (aux - mean(aux(:)))/std(aux(:));
                        oriPinkMatrix(:,colRange) = ApplyOriFilter1(oriFilter(:,:, b)', OFsize(:,:, b), aux);
                    end
                    noiseTex   = Screen('MakeTexture', dpP.window,  prm.noiseSTDmult*noiseMatrix,      [], [], 1);
                    gaborTex   = Screen('MakeTexture', dpP.window,  prm.gaborSTDmult*txP.gabor.matrix, [], [], 1);
                    oriPinkTex = Screen('MakeTexture', dpP.window,  prm.stimSTDmult *oriPinkMatrix,    [], [], 1);

                    xFix = [-crossSize_px/2 crossSize_px/2 0 0];
                    yFix = [0 0 -crossSize_px/2 crossSize_px/2];
                    fixCoords = [xFix; yFix];
                    fixAcquired = false;
                    while ~fixAcquired && keepGoingTrials && ~restartTrial
                        Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        if wasRecording, Eyelink('StopRecording'); end
                        Eyelink('SetOfflineMode');
                        WaitSecs(.1); EyelinkDoDriftCorrection(tkP.el);
                        WaitSecs(.1); Screen('DrawLines', dpP.window, fixCoords, prm.lineWidth_px, drP.white, fixCenters(:, idx, b)', 2);
                        FSonset = Screen('Flip', dpP.window); FPonset = FSonset;
                        Eyelink('StartRecording'); wasRecording = true;
                        Eyelink('Message',sprintf(prm.msg.on.trl{1}, i, nTrials));
                        Eyelink('Message',prm.msg.on.P1);
                        fixCenter = fixCenters(:,idx, b);
                        while true
                            check = false;
                            if (GetSecs - FSonset) > prm.maxCrossDur
                                break;
                            end
                            damn = Eyelink('CheckRecording');
                            if(damn ~= 0), break; end
    
                            if Eyelink('NewFloatSampleAvailable') > 0
                                evt = Eyelink('NewestFloatSample');
                                x_gaze = evt.gx(tkP.Eye);
                                y_gaze = evt.gy(tkP.Eye);
                                check = true;
                            end
                            if check
                                % Se o olho (ou cursor) estiver perto da cruz 
                                % por tempo suficiente, prossegue
                                if vecnorm([x_gaze; y_gaze] - fixCenter) <= minFixDist1
                                    if (GetSecs - FPonset) >= prm.minFixTime1
                                        fixAcquired = true;
                                        break; 
                                    end
                                % Se estiver distante, reinicia a contagem
                                elseif vecnorm([x_gaze; y_gaze] - fixCenter) > minFixDist1
                                    FPonset = GetSecs;
                                end
                            end
                            WaitSecs(.01);
                            [keyIsDown, ~, keyCode] = KbCheck;
                            if keyIsDown
                                if keyCode(escapeKey)
                                    KbReleaseWait;
                                    [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, txP, dpP, drP, prm, 'trial', 0, 2, targetOri(b));
                                    break;
                                end
                            end
                        end
                        if ~fixAcquired && keepGoingTrials && ~restartTrial
                            Eyelink('Message',prm.msg.err.P1);
                            Eyelink('StopRecording');
                            Eyelink('SetOfflineMode');
                            disp('Erro: sem fixação no tempo necessário')
                                
                            tStart = GetSecs;
                            while true
                                alpha = min((GetSecs - tStart) / prm.fadeInDur1, 1);

                                Screen('DrawLines', dpP.window, fixCoords, prm.lineWidth_px, drP.white, fixCenters(:, idx, b)', 2);

                                Screen('FillRect', dpP.window, [overlayColor alpha*.75*drP.white], overlayRect);
                                Screen('TextSize', dpP.window, prm.textSizeEnormous); Screen('TextStyle', dpP.window, 1);
                                DrawFormattedText(dpP.window, 'TEMPO ESGOTADO', 'center', dpP.winRect(4)/2 - 45, [textColor1 alpha*drP.white]);

                                msg = [
                                    'Se o aparelho descalibrou, por favor \n'...
                                    'chame o experimentador.\n\n' ...
                                    'Do contrário, desconsidere a mensagem e \n'...
                                    'aperte ESPAÇO para reiniciar a tentativa'
                                ];
                                Screen('TextSize', dpP.window, prm.textSizeBig); Screen('TextStyle', dpP.window, 0);
                                DrawFormattedText(dpP.window, msg, 'center', dpP.winRect(4)/2 + 15, [textColor1 alpha*drP.white]);
                                Screen('Flip', dpP.window);
                            
                                if alpha >= 1, break; end
                            end
                            WaitSecs(.01)
                            while true
                                [keyIsDown, ~, keyCode] = KbCheck;
                                KbReleaseWait;
                                if keyIsDown
                                    if keyCode(rKey)
                                        Eyelink('Message',prm.msg.pse{5});
                                        EyelinkDoTrackerSetup(tkP.el);
                                        restartTrial = true;
                                        break;
                                    elseif keyCode(spaceKey)
                                        disp('Prossegue sem recalibragem')
                                        break;
                                    elseif keyCode(escapeKey)
                                        [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, txP, dpP, drP, prm, 'trial', 0, 2, targetOri(b));
                                        break;
                                    end
                                end
                                WaitSecs(0.02);
                            end
                        end
                    end
                    
%% Fase 2: tela de estímulos
                    blinkIdx = setdiff(1:nStims, fixIdx(b, idx));
                    vbl = Screen('Flip', dpP.window);
                    if keepGoingTrials && ~restartTrial
                        alphas = ones(nStims, 1);
                        foragingDrawMain(dpP.window, gaborTex, noiseTex, srcRects, dstRects, orientation(:, idx, b), txP, [repmat(alphas', [3,1]); ones(1, nStims)]);
                        P2On = Screen('Flip', dpP.window, vbl + 0.5 * dpP.ifi);
                        Eyelink('Message',prm.msg.on.P2);
                    end

    
%% Fase 3: tela com ruído rosa
                    if ~restartTrial && keepGoingTrials
                        Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
                        Screen('DrawTextures', dpP.window, oriPinkTex, srcRects(:,blinkIdx), dstRects(:,blinkIdx), orientation(blinkIdx, idx, b), [], [], [], []);

                        Screen('DrawTextures', dpP.window, gaborTex, [], dstRects(:,fixIdx(b, idx)), orientation(fixIdx(b, idx), idx, b));
                        Screen('DrawTextures', dpP.window, noiseTex, srcRects(:,fixIdx(b, idx)), dstRects(:,fixIdx(b, idx)), orientation(fixIdx(b, idx), idx, b));
                        
                        Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                        Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects(:, fixIdx(b, idx)), orientation(fixIdx(b, idx), idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
                        
                        Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects(:, blinkIdx), orientation(blinkIdx, idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
                        Screen('Close', oriPinkTex); Screen('Close', noiseTex); Screen('Close', gaborTex);
                        Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                        P3On = Screen('Flip', dpP.window, P2On + jitterTimes(b, idx));
                        
                        Eyelink('Message',prm.msg.on.P3);

                        tNow = GetSecs;
                        keptFixP3 = false; 
                        tic
                        while tNow - P3On < prm.pinkNoiseDur
                            check = false;
                            damn = Eyelink('CheckRecording');
                            if(damn ~= 0), break; end
                            if Eyelink('NewFloatSampleAvailable') > 0
                                evt = Eyelink('NewestFloatSample');
                                x_gaze = evt.gx(tkP.Eye);
                                y_gaze = evt.gy(tkP.Eye);
                                check = true;
                            end
        
                            if check
                                keptFixP3 = vecnorm([x_gaze; y_gaze] - fixCenters(:, idx, b)) <= minFixDist3;
                                if ~keptFixP3, break; end
                            end
                            WaitSecs(.0005);
                            tNow = GetSecs;
                        end
                        toc

                        Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                        PMOn = Screen('Flip', dpP.window);
                        Eyelink('Message',prm.msg.off.P3);
                        
%% Fase PM: pedestais
                        Eyelink('Message',prm.msg.on.PM);
                        keptFixPM = false;  
                        tNow = GetSecs;
                        while tNow - PMOn < prm.postModDurStair
                            check = false;
                            damn = Eyelink('CheckRecording');
                            if(damn ~= 0), break; end
                            if Eyelink('NewFloatSampleAvailable') > 0
                                evt = Eyelink('NewestFloatSample');
                                x_gaze = evt.gx(tkP.Eye);
                                y_gaze = evt.gy(tkP.Eye);
                                check = true;
                            end
        
                            if check
                                keptFixPM = vecnorm([x_gaze; y_gaze] - fixCenters(:, idx, b)) <= minFixDist3;
                                if ~keptFixPM, break; end
                            end
                            WaitSecs(.0005);
                            tNow = GetSecs;
                        end

                        WaitSecs(prm.fadeInDelay1);

                        Eyelink('Message',prm.msg.off.PM);
                        
%% Fase 4: tela de reporte
                        allTargets = nan(1,nStims); allColors2 = drP.allColors;
                        Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                        % Os demais pontos que não o fixado antes da atual
                        nbhd = NeighborsOrder(stimCenters(:, :, idx, b), fixIdx(b, idx));

                        stimsToReport = nbhd(1:3);
                        orderToReportStims = stimsToReport(randperm(3));
                        Eyelink('Message',prm.msg.on.P4);
                        
                        for j=1:length(orderToReportStims)
                            KbReleaseWait;
                            rectColors = allColors2; rectColors(:, orderToReportStims(j)) = drP.orange;
                        
                            rectPW = drP.allPW; rectPW(:, orderToReportStims(j)) = prm.pW2;
                            currTarget = -1;
                            abort  = false;
                        
                            while ~abort
                                [keyIsDown, ~, keyCode] = KbCheck;
                                if keyIsDown
                                    if keyCode(nonTargetKey)
                                        currTarget = 0;
                                    elseif keyCode(targetKey)
                                        currTarget = +1;
                                    elseif keyCode(spaceKey)
                                        KbReleaseWait;
                                        if currTarget ~= -1
                                            abort = true;
                                            allColors2(:, orderToReportStims(j)) = drP.whiteGrey;
                                        end
                                    end
                                end
                                allTargets(orderToReportStims(j)) = currTarget;
                                Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                                foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, rectColors, allTargets, targetOri(b), rectPW);
                            end
                        end
                        
                        Eyelink('Message',prm.msg.off.P4);

                        feedback = (orientation(:,idx,b) == targetOri(b))' + allTargets;
                        feedback(rem(feedback,2) == 0) = 2; feedback(rem(feedback,2) == 1) = 0;
                        feedback = feedback/2;

                        trialFeedback{b, i} = [orderToReportStims; zeros(1, 3); feedback(orderToReportStims)];
                        
                        Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                        % foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW);
                        foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW, feedback, drP.red, drP.green);
                        WaitSecs(.5);

                        trialOrder(2, i, b) = 1;

                        i = i + 1;
                        trialIdxUp = true;

                        restartTrial = ~(keptFixP3 && keptFixPM);
                    end

                    if ~keepGoingTrials || restartTrial
                        fprintf('Trial ruim, pois ~keepGoingTrials = %d e restartTrial = %d\n keptFixP3 = %d, keptFixPM = %d\n', ~keepGoingTrials, restartTrial, keptFixP3, keptFixPM)
                        Eyelink('Message',prm.msg.err.trl{1});
                        retryCount(trialQueue(i)) = retryCount(trialQueue(i)) + 1;
                        if keepGoingTrials
                            if retryCount(trialQueue(i)) > prm.maxRetries
                                Eyelink('Message',prm.msg.err.trl{2});
                                warning('Trial %d excede o máximo de repetições. Prosseguindo', trialQueue(i));
                                i = i + 1;
                                trialIdxUp = true;
                            else
                                trialQueue(i:end) = [trialQueue(i+1:end) trialQueue(i)];
                            end
                        end
                    else
                        disp('Trial bom')
%% Atualiza o staircase se o trial foi útil
                        for j=1:length(orderToReportStims)
                            RF(b) = PAL_AMRF_updateRF(RF(b), -aSigma(b), feedback(orderToReportStims(j)));
                        end
                        aSigma(b) = -RF(b).mean;
                        fprintf('Atualizo o RF do bloco %d, apresentamos: ', b); disp(RF(b).x)
                        fprintf('\nA média evoluiu como: '); disp(RF(b).xStaircase)
                        fprintf('\nPor isso o novo aSigma é %.4f\n', aSigma(b));

                        [auxOriFilter, auxOFsize] = MakeOriFilter1(txP.gabor.size_px, aSigma(b), prm.rSigma2);
                        oriFilter(:,:,b) = auxOriFilter;
                        OFsize(:,:,b)    = auxOFsize;
                    end
                    trialOffset = GetSecs; %#ok<NASGU>
                    Eyelink('Message', sprintf(prm.msg.off.trl{1}, i - trialIdxUp, nTrials));
                    Screen('Flip', dpP.window);
                    WaitSecs(0.1);
                    Eyelink('SetOfflineMode');
                    Eyelink('StopRecording');
                end
                if restartBlock
                    Eyelink('Message',prm.msg.err.blk);
                    Eyelink('Message',sprintf(prm.msg.off.blk{1}, b, nBlocks));
                    restartBlock = false;
                else
                    blockOffset = GetSecs; %#ok<NASGU>
                    Eyelink('Message',sprintf(prm.msg.off.blk{1}, b, nBlocks));
                    b = b+1;
                end
            end
            Eyelink('Message',prm.msg.off.stc);
            resultsStair.fixCenters = fixCenters;
            resultsStair.stimCenters = stimCenters;
            resultsStair.orientation = orientation;
            resultsStair.nTs = nTs;
            resultsStair.targetOri = targetOri;
            resultsStair.nStimsToReport = nStimsToReport;
            resultsStair.trialOrder = trialOrder;
            resultsStair.trialFeedback = trialFeedback;
            resultsStair.aSigma    = aSigma;
            resultsStair.oriFilter = oriFilter;
            resultsStair.OFsize   = OFsize;
            resultsStair.staircase   = RF;
        catch
            if ~exist('fixCenters', 'var'),   fixCenters = []; end
            if ~exist('stimCenters', 'var'),  stimCenters = []; end
            if ~exist('orientation', 'var'),  orientation = []; end
            if ~exist('nTs', 'var'),       nTs = []; end
            if ~exist('targetOri', 'var'), targetOri = []; end
            if ~exist('nStimsToReport', 'var'),     nStimsToReport = []; end
            if ~exist('orderToReportSets', 'var'),  orderToReportSets = []; end
            if ~exist('trialOrder', 'var'),      trialOrder = []; end
            if ~exist('trialFeedback', 'var'),   trialFeedback = []; end

            resultsStair.fixCenters = fixCenters;
            resultsStair.stimCenters = stimCenters;
            resultsStair.orientation = orientation;
            resultsStair.nTs = nTs;
            resultsStair.targetOri = targetOri;
            resultsStair.nStimsToReport = nStimsToReport;
            resultsStair.orderToReportSets = orderToReportSets;
            resultsStair.trialOrder = trialOrder;
            resultsStair.trialFeedback = trialFeedback;
            
%             foragingSave(tkS, 2, prm, dpP, drP, tkP, txP, results);
            cleanup(dpP.window);
            diary off;
            psychrethrow(psychlasterror);
        end
        if b == nBlocks && keepGoingBlocks
            plotStaircase(RF)
        end
end
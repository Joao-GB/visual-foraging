function [resultsStair, tkS] = runStaircase(tkP, dpP, drP, txP, prm, mode, tkS)
% Nessa versão, usa um único staircase para todos os tipos de estímulo, 
% alimentando o ajuste da curva
        Screen('Flip', dpP.window);
        nBlocks = prm.nBlocksStair; nTrials = prm.nTrialsStair;

        if isfield(tkP,'pinkNoiseDur'), prm.pinkNoiseDur = tkP.pinkNoiseDur; end

        if mode <= 3
            prm.pinkNoiseDur = prm.cursorPinkNoiseDur;
            nTrials = prm.nTrialsStairTrain;
            prm.aSigma = prm.aSigmaTrain;
            tkP.stairPrev = [];
        
        end
        
        leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4}; rKey = tkP.keys{5};
        if strcmp(tkP.targetKey, 'right'), targetKey = rightKey; nonTargetKey = leftKey;
        else,                              targetKey = leftKey; nonTargetKey = rightKey;
        end
        % Não me interessam os outputs de distribuições temporais, pois no
        % staircase o tempo não será por quantidade de fixações
        nTrialsBuffered = nTrials + prm.nBufferTrials;
        nStims = prm.nStimsStair;
        [nTs,  ~, targetOri, ~, ~, ~] = getForagingDistributions1(nStims, tkP.nMinFix, tkP.nMaxFix, nTrialsBuffered, nBlocks, prm);
        
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

        auxidx = zeros(1, nBlocks);
        stimRange = (-prm.sigmaMax:prm.sigmaStep:-prm.sigmaMin);
        if isscalar(prm.aSigma), prm.aSigma = repmat(prm.aSigma, 1, nBlocks); end
        for b=1:nBlocks
            for i=1:nTrialsBuffered
                [currFixCenter, currStimCenter, ~, fI] = getStimLocations2_1(dpP.winRect(3:4), [dpP.winCenter 1], tkP.nStims, minDist_px, txP.gabor.size_px);
                fixCenters(:, i, b)     = currFixCenter;
                aux = NeighborsOrder(currStimCenter, fI, []);
                stimCenters(:, :, i, b) = currStimCenter(:, [aux(1:nStims-1), fI]);
                fixIdx(b, i) = nStims;
                orientation(randperm(nStims, nTs(b, i)), i, b) = targetOri(b);
            end
            [~, auxidx(b)] = min(abs(stimRange - (-prm.aSigma(prm.allOriMap(targetOri(b))))));
            prm.aSigma(b) = stimRange(auxidx(b));
        end
    
        minFixDist3 = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.fixROIradius3_dva);

%% Cria o objeto que registra todo o staircase

        
        if tkP.stairBurnIn
            burninTrials = prm.burninTrials;
            auxStimRange = -flip(stimRange);
            burnInASigma = linspace(prm.sigmaMin, prm.sigmaMax, burninTrials);
            burnInASigma = interp1(auxStimRange, auxStimRange, burnInASigma, 'nearest');
            burnInASigma = burnInASigma(randperm(burninTrials));
            aSigma = ones(1, nBlocks)*burnInASigma(1);
        else
            burninTrials = 0;          % Neutraliza o loop posterior
            burnInASigma = prm.aSigma; aSigma = prm.aSigma;
        end
        PF = @PAL_CumulativeNormal;

        % Como apenas gama é determinado, definimos os demais parâmetros
        % através de intervalos
        gamma = 0.5;
        marginalize = 4;
        
        priorAlphaRange = linspace(-prm.sigmaMax, -prm.sigmaMin, prm.grainAlpha);
        priorBetaRange  = 10.^(linspace(log10(prm.betaMin), log10(prm.betaMax), prm.grainBeta));
        priorLambdaRange = (0:0.01:0.1) + .001;
            
        % Veja que removi a condicional que usa burn-in
        for b = 1:nBlocks
            PM(b) = PAL_AMPM_setupPM(...
            'priorAlphaRange', single(priorAlphaRange), ...
            'priorBetaRange', single(priorBetaRange), ...
            'priorGammaRange', single(gamma), ...
            'priorLambdaRange', single(priorLambdaRange), ...
            'numtrials', nBlocks * nTrials * (nStims-1), ...
            'PF', PF, ...
            'stimRange', single(stimRange), ...
            'marginalize', marginalize);

            % Construção do prior customizado usando matrizes multidimensionais do PM
            if tkP.stairBurnIn == 0 && isfield(tkP, 'stairPrev') && ~isempty(tkP.stairPrev)
                k = find([tkP.stairPrev.tgtOri] == targetOri(b),1);
                prevSessionAlpha = tkP.stairPrev(k).threshold(end);
                prevSessionBeta  = tkP.stairPrev(k).slope(end);

                prior = PAL_pdfNormal(PM(b).priorAlphas, prevSessionAlpha, prm.priorStdStair2);
                prior = prior .* PAL_pdfNormal(PM(b).priorBetas, prevSessionBeta, prm.priorBetaStdStair2);
            else
                prior = PAL_pdfNormal(PM(b).priorAlphas, -prm.priorMeanStair, prm.priorStdStair);
                prior = prior .* PAL_pdfNormal(PM(b).priorBetas, prm.priorBetaMeanStair, prm.priorBetaStdStair);
            end

            prior = prior .* PAL_pdfBeta(PM(b).priorLambdas, prm.priorLambdaMeanStair, prm.priorLambdaStdStair, 'meanandconcentration');

            % Normaliza a grade de probabilidade
            prior = prior ./ sum(prior(:));
            PM(b) = PAL_AMPM_setupPM(PM(b), 'prior', prior); %#ok<*AGROW>
        end
        for b = 1:nBlocks
            PM(b).tgtOri = targetOri(b);
        end
        startPM = PM;
        startASigma = aSigma;

        [oriFilter, OFsize] = MakeOriFilter1(txP.gabor.size_px, burnInASigma(1), prm.rSigma2);

        oriFilter = repmat(oriFilter, [1, 1, nBlocks]);
        OFsize    = repmat(OFsize, [1, 1, nBlocks]);
        jitterTimes = rand(nBlocks, nTrialsBuffered)*(prm.maxJitterStair-prm.minJitterStair)+prm.minJitterStair;

        suspend = 0;

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
            
            stimHistory = cell(1,nBlocks);
            respHistory = cell(1,nBlocks);

            while b <= nBlocks && keepGoingBlocks
%                 if b > 1
%                     startPM(b) = PM(b-1);
%                     startASigma(b) = aSigma(b-1);
%                 end
                stimHistory{b} = [];
                respHistory{b} = [];
                PM(b) = startPM(b);
                aSigma(prm.allOriMap(targetOri(b))) = startASigma(prm.allOriMap(targetOri(b)));
                [~, auxidx(b)] = min(abs(stimRange - (-startASigma(prm.allOriMap(targetOri(b))))));

                [auxOriFilter, auxOFsize] = MakeOriFilter1(txP.gabor.size_px, aSigma(prm.allOriMap(targetOri(b))), prm.rSigma2);
                oriFilter(:,:,b) = auxOriFilter;
                OFsize(:,:,b)    = auxOFsize;

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
                    keptFixP3 = false;
                    keptFixPM = false;
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
                        nbhd = NeighborsOrder(stimCenters(:, :, idx, b), fixIdx(b, idx), []);

                        stimsToReport = nbhd(1:nStims-1);
                        orderToReportStims = stimsToReport(randperm(nStims-1));
                        WaitSecs(.1);
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

                        trialFeedback{b, i} = [orderToReportStims; zeros(1, nStims-1); feedback(orderToReportStims)];

                        restartTrial = ~(keptFixP3 && keptFixPM);
                        
                        Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                        


                        if restartTrial
                            warningFlip(dpP.window, stimCenters(:, :, idx, b), resizeRect(dstRects, .25), fixIdx(b, idx), txP.gabor.size_px, drP.allPW, drP.darkRed);
                        else
                            trialOrder(2, i, b) = 1;
                            i = i + 1;
                            trialIdxUp = true;

                            if mode <= 3
                                foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW, feedback, drP.red, drP.green);
                            else
                                foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW);
                            end
                        end

                        WaitSecs(.5);
                        
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
                        if aSigma(prm.allOriMap(targetOri(b))) <= (prm.sigmaMin + .5*prm.sigmaRem) && prm.avoidConsecutive, suspend = 1;
                        end
                        if suspend == 1, suspend = rand(1) > 1./prm.stairWaitTime;
                        end
                        
                        for j=1:length(orderToReportStims)
                            PM(b) = PAL_AMPM_updatePM(PM(b), feedback(orderToReportStims(j)), 'xIndex', auxidx(b), 'fixLapse', suspend);
                            stimHistory{b} = [stimHistory{b} stimRange(auxidx(b))];
                            respHistory{b} = [respHistory{b} feedback(orderToReportStims(j))];
                        end

                        fprintf('Atualizo o PM do bloco %d, apresentamos: %.4f (no xCurrent: %.4f)\n', b, aSigma(prm.allOriMap(targetOri(b))), PM(b).xCurrent);
                        % fprintf('\nA média evoluiu como: '); disp(PM(b).xStaircase)

                        % Quando chega nessa parte, o i já foi incrementado 
                        % na parte com trialIdxUp. Por iso i em vez de i+1
                        if i <= burninTrials && b == 1
                            [~, auxidx(b)] = min(abs(stimRange - (-burnInASigma(i))));
                            aSigma(prm.allOriMap(targetOri(b))) = stimRange(auxidx(b));
                            fprintf('\nMas como é burn-in, o novo aSigma é %.4f\n', aSigma(prm.allOriMap(targetOri(b))));
                            % PM(b).xCurrent = -aSigma(prm.allOriMap(targetOri(b)));
                        else
                            [~, auxidx(b)] = min(abs(stimRange - (PM(b).xCurrent)));
                            aSigma(prm.allOriMap(targetOri(b))) = -PM(b).xCurrent;
                            fprintf('\nPor isso o novo aSigma é %.4f\n', aSigma(prm.allOriMap(targetOri(b))));
                        end

                        currentAlphaEst = PM(b).threshold(end);
                        currentBetaEst  = PM(b).slope(end);
                        fprintf('Limiar atual estimado (Alpha): %.2f\n', currentAlphaEst);
                        fprintf('Inclinação atual estimada (Beta): %.2f\n\n', currentBetaEst);

                        [auxOriFilter, auxOFsize] = MakeOriFilter1(txP.gabor.size_px, aSigma(prm.allOriMap(targetOri(b))), prm.rSigma2);
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
            resultsStair.trialOrder = trialOrder;
            resultsStair.trialFeedback = trialFeedback;

            resultsStair.aSigma75  = aSigma;

            newLevelASigma = ones(1, nBlocks)*prm.aSigma(1);

            if b == nBlocks + 1 && keepGoingBlocks
            
                searchGrid.alpha  = priorAlphaRange;
                searchGrid.beta   = priorBetaRange;
                searchGrid.gamma  = gamma;
                searchGrid.lambda = priorLambdaRange;
                
                paramsFree = [1 1 0 1]; 
                for i = 1:nBlocks
                    % % MLE restrito ao grid
                    % [~, maxIndex] = PAL_findMax(PM(i).pdf);
                    % MLalpha = priorAlphaRange(maxIndex(1));
                    % MLbeta = priorBetaRange(maxIndex(2));
                    % MLlambda = priorLambdaRange(maxIndex(4));
    
    
                    % MLE preciso
                    fprintf('Comprimentos dos históricos de estímulos e respostas para bloco %d: %d e %d\n', i, numel(stimHistory{i}), numel(respHistory{i}));
                    StimLevels = unique(stimHistory{i});
    
                    NumPos   = zeros(size(StimLevels));
                    OutOfNum = zeros(size(StimLevels));
                    
                    for k = 1:numel(StimLevels)
                        idx = stimHistory{i} == StimLevels(k);
                        OutOfNum(k) = sum(idx);
                        NumPos(k)   = sum(respHistory{i}(idx));
                    end
    
                    [paramsValues, ~] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, paramsFree, PF);
                    MLalpha  = paramsValues(1);
                    MLbeta   = paramsValues(2);
                    MLlambda = paramsValues(4);
    
    
                    fprintf('ANTES\nAlpha: ML = %.3f e final = %.3f\n', MLalpha, PM(i).threshold(end));
                    fprintf('Beta: ML = %.3f e final = %.3f\n', MLbeta, PM(i).slope(end));
                    fprintf('Lambda: ML = %.3f e final = %.3f\n\n', MLlambda, PM(i).lapse(end));
    
                    if isnan(MLalpha) || isempty(MLalpha), MLalpha = PM(i).threshold(end); end
                    if isnan(MLbeta) || isempty(MLbeta), MLbeta  = PM(i).slope(end); end
                    if isnan(MLlambda) || isempty(MLlambda), MLlambda  = PM(i).lapse(end); end

                    fprintf('DEPOIS\nAlpha: ML = %.3f e final = %.3f\n', MLalpha, PM(i).threshold(end));
                    fprintf('Beta: ML = %.3f e final = %.3f\n', MLbeta, PM(i).slope(end));
                    fprintf('Lambda: ML = %.3f e final = %.3f\n\n', MLlambda, PM(i).lapse(end));
    
                    PM(i).threshold(end) = MLalpha;
                    % PM(i).slope(end)     = MLbeta;
                    PM(i).lapse(end)     = MLlambda;
                    
                    if mode <= 3 || isempty(PM(i).threshold)
                        newLevelASigma(prm.allOriMap(targetOri(i))) = aSigma(prm.allOriMap(targetOri(i)));
                    else
                        newLevelASigma(prm.allOriMap(targetOri(i))) = -PAL_CumulativeNormal([MLalpha, MLbeta, gamma, MLlambda], min(prm.stairLevel, 1-MLlambda-.0001), 'inverse');
                    end
                end
            end

            resultsStair.aSigma    = newLevelASigma;
            resultsStair.oriFilter = oriFilter;
            resultsStair.OFsize    = OFsize;
            resultsStair.staircase = rmfield(PM, {'priorAlphas', 'priorBetas', 'priorGammas', 'priorLambdas', 'priorModels', 'LUT', 'posteriorTplus1givenSuccess', 'posteriorTplus1givenFailure'});
        catch
            if ~exist('fixCenters', 'var'),   fixCenters = []; end
            if ~exist('stimCenters', 'var'),  stimCenters = []; end
            if ~exist('orientation', 'var'),  orientation = []; end
            if ~exist('nTs', 'var'),       nTs = []; end
            if ~exist('targetOri', 'var'), targetOri = []; end
            if ~exist('trialOrder', 'var'),      trialOrder = []; end
            if ~exist('trialFeedback', 'var'),   trialFeedback = []; end
            if ~exist('aSigma', 'var'), aSigma = []; end
            if ~exist('newLevelASigma', 'var'),      newLevelASigma = []; end
            if ~exist('oriFilter', 'var'),   oriFilter = []; end
            if ~exist('OFsize', 'var'),      OFsize = []; end
            if ~exist('PM', 'var'),   PM = [];
            else
                PM = rmfield(PM, {'priorAlphas', 'priorBetas', 'priorGammas', 'priorLambdas', 'priorModels', 'LUT', 'posteriorTplus1givenSuccess', 'posteriorTplus1givenFailure'});
            end

            resultsStair.fixCenters = fixCenters;
            resultsStair.stimCenters = stimCenters;
            resultsStair.orientation = orientation;
            resultsStair.nTs = nTs;
            resultsStair.targetOri = targetOri;
            resultsStair.trialOrder = trialOrder;
            resultsStair.trialFeedback = trialFeedback;
            resultsStair.aSigma75 = aSigma;
            resultsStair.aSigma = newLevelASigma;
            resultsStair.oriFilter = oriFilter;
            resultsStair.OFsize = OFsize;
            resultsStair.staircase = PM;
            
            cleanup(dpP.window);
            diary off;
            psychrethrow(psychlasterror);
        end
        clearvars -except b resultsStair nBlocks keepGoingBlocks mode tkP dpP drP prm PM aSigma newLevelASigma targetOri tkS
        if b == nBlocks + 1 && keepGoingBlocks && mode > 3
            inspectStaircase(tkP, dpP, drP, prm, PM, aSigma, newLevelASigma, targetOri);
            tkS(1,2) = 1;
        end
end
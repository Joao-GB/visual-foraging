function [tkP, tkS, results] = trainPreSaccForaging(tkP, dpP, drP, txP, prm, debug, mode, tkS)
        if isempty(mode), mode = 5; end

        Screen('Flip', dpP.window);

        if isfield(tkP,'pinkNoiseDur'), prm.pinkNoiseDur = tkP.pinkNoiseDur; end

        suffix = prm.msg.suffix{mode};
        
        leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4}; rKey = tkP.keys{5};
        if strcmp(tkP.targetKey, 'right')
            targetKey = rightKey; nonTargetKey = leftKey;
        else
            targetKey = leftKey; nonTargetKey = rightKey;
        end

        % Inicializa o emaFix com a mediana da fila de fixações
        if ~isfield(tkP.fixProps, 'emaFix')
            tkP.fixProps.emaFix = median(tkP.fixQueue);
        end
    
    %% 1)

        nBlocks = tkP.nBlocks;
        nTrials = tkP.nTrials;
        L = numel(prm.allOri);
        tkP.nBlocks = L;
        tkP.nTrials = prm.nTrialsTrain3;

        fprintf('Tempo inicial: %.5f\n', toc)
        tic;
    %% 2) Define as distribuições das condições dos trials e dos blocos
        nTrialsBuffered = tkP.nTrials + prm.nBufferTrials;
        tkP.nStims = 3;
        [nTs, nStims, ~, modTimes, nStimsToReport, orderToReportSets] = getForagingDistributions1(tkP.nStims, tkP.nMinFix, tkP.nMaxFix, nTrialsBuffered, tkP.nBlocks, prm);
        targetOri = prm.allOri(randperm(tkP.nBlocks));
        modTimes(:) = 1;

        drP.allColors = drP.white*ones(3, nStims);
        drP.allPW     = prm.pW1*ones(1,nStims);
        
        fprintf('Tempo getDistr: %.5f\n', toc)
        tic;
    %% 3) Define matrizes usadas para os estímulos
        % (a) Matrizes com orientação e centros dos estímulos. Sem os alvos,
        %     há distratores com orientações aleatórias
        %     (cf. (e) para adição de alvos)
        stimCenters = zeros(2, nStims, nTrialsBuffered, tkP.nBlocks);
        orientation = zeros(nStims, nTrialsBuffered, tkP.nBlocks);
        for b=1:tkP.nBlocks
            lowerBound = targetOri(b) + prm.nbhdRadius;
            upperBound = 180 + (targetOri(b) - prm.nbhdRadius);
            orientation(:,:,b) = mod(rand(nStims, nTrialsBuffered)*(upperBound-lowerBound)+lowerBound, 180);
        end
        fprintf('Tempo orientações distratoras: %.5f\n', toc)
        tic;
        % (b) Tamanho e matrizes para as cruzes de fixação 
        crossSize_px = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.crossSize_dva);
        fixCenters  = zeros(2, nTrialsBuffered, tkP.nBlocks);

        % (c) Matrizes com os ruído e o centro dos retângulos a serem
        %     plotados (srcRect)
        noiseMatrix = zeros(txP.gabor.size_px, (nStims)*txP.gabor.size_px);
        oriPinkMatrix = zeros(txP.gabor.size_px, (nStims)*txP.gabor.size_px);

        fprintf('Tempo matrizes nulas: %.5f\n', toc)
        tic;
        noiseCenters = [0:txP.gabor.size_px:(nStims-1)*txP.gabor.size_px; zeros(1, nStims)] + txP.gabor.size_px/2;
        baseRect = [0 0 txP.gabor.size_px txP.gabor.size_px];

        srcRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', noiseCenters(1,:), noiseCenters(2,:));

        % Já que gaborSize_dva é o diâmetro do Gabor, os centros precisam
        % distar pelo menos 1 diâmetro mais a minDist
        minDist_px   = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.minDist_dva+prm.gaborSize_dva);
        minFixDist1 = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.fixROIradius1_dva);

        fprintf('Tempo rects a centers: %.5f\n', toc);
        tic;
        for b=1:tkP.nBlocks
            for i=1:nTrialsBuffered
        % (d) Matrizes com os centros dos estímulos (i.e., centros dos dstRects)
        %     e das cruzes de fixação
                [currFixCenter, currStimCenter, ~] = getStimLocations2_1(dpP.winRect(3:4), [dpP.winCenter 1], nStims, minDist_px, txP.gabor.size_px);
                fixCenters(:, i, b)     = currFixCenter;
                stimCenters(:, :, i, b) = currStimCenter;
        % (e) Matriz com orientações tem os alvos adicionados
                orientation(randperm(nStims, nTs(b, i)), i, b) = targetOri(b);
            end
        end

        % (f) Distância mínima para considerar fixação em alvo
        % pós-modificação (raio, não diâmetro)
         minFixDist3 = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.fixROIradius3_dva);

        auxFixQueue = zeros(1, nStims);
        fprintf('Tempo preencher locais de sestímulos e orientações de alvos: %.5f\n', toc);
    %% 4) Início dos blocos e trials
        fprintf('----Início da sessão (pré-sacádica) ----\n')
        try
            Screen('TextFont', dpP.window, prm.textFont);
            keepGoingBlocks = true;
            wasRecording = false;
            restartBlock = false;
            blocksCompleted = false;
            b = 1;

            trialOrder = zeros(2, tkP.nTrials, tkP.nBlocks);

            % O vetor guarda os índices dos estímulos em destaque na fase 4 (i.e., sobre os quais o sujeito tinha que responder)
            % na ordem em que foram perguntados, identificando em 3 linhas
            %   i. o índice do estímulo;
            %  ii. se era alvo (0), já visto (-1) ou ainda não (1)
            % iii. se acertou ou errou
            trialFeedback = cell(tkP.nBlocks, tkP.nTrials);
            orderToReportMap = [-1 0 1];

            seenStimsQueue = cell(tkP.nBlocks, tkP.nTrials);

            nSnbhd = zeros(tkP.nBlocks, tkP.nTrials);

            isSaccSeen    = nan(tkP.nBlocks, tkP.nTrials);
            isP3earlyStop = nan(tkP.nBlocks, tkP.nTrials);

            Eyelink('Message',sprintf(prm.msg.on.ses{1}, suffix));
            EyelinkDoTrackerSetup(tkP.el);

            while b <= tkP.nBlocks && keepGoingBlocks
                trialOrder(:,:,b) = 0;
                
                % Para as telas de erro
                overlayRect = CenterRectOnPointd([dpP.winRect(1:3) dpP.winRect(4)/2.5], dpP.winRect(3)/2, dpP.winRect(4)/2);
                overlayColor = repmat(drP.whiteGrey, 1, 3);
                textColor1   = repmat(drP.black, 1, 3);
                textColor2   = repmat(.5*(drP.grey+drP.black), 1, 3);

        % (a) Mensagem com orientação a ser procurada
                repeatMessage = true;
                
                while repeatMessage && keepGoingBlocks
                    Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
                    Screen('DrawTexture', dpP.window, txP.exampleGabor.tex, [], [], targetOri(b), [], [], [], [], []);
                    Screen('DrawTexture', dpP.window, txP.exampleNoise.tex, [], [], [], [], [], [], [], []);
                    Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                    Screen('DrawTextures', dpP.window, txP.exampleBlob.tex, [], [], targetOri(b), [], [], [0 0 0 1]', [], [], txP.exampleBlob.props); 
        
                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    Screen('TextSize', dpP.window, prm.textSizeBigger);
                    countBlockText  = sprintf('Bloco %d de %d', b, tkP.nBlocks);
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
                        [keepGoingBlocks, ~, ~, ~] = pauseHandle(keepGoingBlocks, [], [], [], wasRecording, tkP, txP, dpP, drP, prm, 'block', debug, mode, []);
                    end
                end

                blockOnset = Screen('Flip', dpP.window); %#ok<NASGU>
                
                if debug == 0 && keepGoingBlocks
                    Eyelink('Message', sprintf(prm.msg.on.blk{1}, b, tkP.nBlocks));
                end
                i = 1;
                % Reordeno os trials para que, caso reinicie o bloco, a
                % ordem seja diferente
                trialQueue = randperm(nTrialsBuffered);
                retryCount = zeros(1, nTrialsBuffered);

                keepGoingTrials = keepGoingBlocks;
                % Esse loop-mor é o único que não considera restartTrial;
                % não quero interroper todos os trials, apenas o atual, e
                % recomeçar este loop
                while i <= tkP.nTrials && keepGoingTrials
                    trialIdxUp = false;
                    seenStimsQueue{b, i} = [];
                    idx = trialQueue(i);
                    trialOrder(1, i, b) = idx;
                    restartTrial = false;
                    fprintf('\nIdx Bloco: %d\n# Trial: %d\n Idx Trial: %d\n', b, i, idx)
                    fprintf('ATENÇÃO: %d/%d visitas até modificar\n', modTimes(b, idx), nStims)

        % (b) Obtém a mediana do vetor de tempos de fixação
                    med  = median(tkP.fixQueue);
                    medFixTime = median(tkP.fixQueue);
                    maxTrialDur = prm.maxTrialDurFactor*tkP.nStims*medFixTime;
        % (c) Cria os retângulos de destino com base nas coordenadas dos
        %     centros dos estímulos
                    dstRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', stimCenters(1,:, idx, b), stimCenters(2,:, idx, b));
                    
                
                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    if debug == 0 && mode >= 2
    
        % iii. Desenha, na tela do Host PC, linhas indicando a orientação de cada estímulo
                        % Preciso do sinal de menos no dy pois as coord.
                        % do Eyelink são diferentes das do PTB
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
                    end
                    
    %% 5) Início Fase 1: tela de fixação
                    % Se estiver no modo de cursor, cria uma textura em janela
                    % offscreen sobre a qual será desenhado o símbolo do cursor
                    if mode == 1
                        bg = Screen('OpenOffscreenWindow', dpP.window, drP.grey, [], 64);
                        auxWin  = bg;
                    else
                        auxWin = dpP.window;
                    end

        % (d) Cria as texturas de gabores e de ruído com e sem orientação
        %     (pois não faz sentido armazená-las), sem desenhá-las

%                     auxNoiseMatrix = butterFilter(pinkNoise(txP.gabor.size_px, (nStims)*txP.gabor.size_px), txP.noiseLoCutFreq, txP.noiseHiCutFreq);
                    auxNoiseMatrix = pinkNoise(txP.gabor.size_px, (nStims)*txP.gabor.size_px);
                    fprintf('Vai usar filtro %.2f para orientaçao %d\n', prm.aSigma(prm.allOriMap(targetOri(b))), targetOri(b));
                    for j=1:(nStims)
                        colRange = ((j-1)*txP.gabor.size_px+1):(j*txP.gabor.size_px);
                        aux = auxNoiseMatrix(:,colRange);
                        aux1 = butterFilter(aux, txP.noiseLoCutFreq, txP.noiseHiCutFreq);
                        noiseMatrix(:, colRange) = (aux1 - mean(aux1(:)))/std(aux1(:));
                        oriPinkMatrix(:,colRange) = ApplyOriFilter1(txP.oriFilter{prm.allOriMap(targetOri(b))}', txP.OFsize{prm.allOriMap(targetOri(b))}, aux);
                    end
                    noiseTex   = Screen('MakeTexture', dpP.window,  prm.noiseSTDmult*noiseMatrix,      [], [], 1);
                    gaborTex   = Screen('MakeTexture', dpP.window,  prm.gaborSTDmult*txP.gabor.matrix, [], [], 1);
                    oriPinkTex = Screen('MakeTexture', dpP.window,  prm.stimSTDmult *oriPinkMatrix,    [], [], 1);

        % (f) Desenha a cruz de fixação
                    xFix = [-crossSize_px/2 crossSize_px/2 0 0];
                    yFix = [0 0 -crossSize_px/2 crossSize_px/2];
                    fixCoords = [xFix; yFix];
                    fixAcquired = false;
                    while ~fixAcquired && keepGoingTrials && ~restartTrial

                        Screen('BlendFunction', auxWin, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        if debug == 0 && mode >= 2
        % i. Deixa de registrar até StartRecording
                            if wasRecording, Eyelink('StopRecording'); end
                            Eyelink('SetOfflineMode');
                            WaitSecs(.1);
        % ii. Faz drift correction
                            EyelinkDoDriftCorrection(tkP.el);
                            WaitSecs(.1);
                        end
                        Screen('DrawLines', auxWin, fixCoords, prm.lineWidth_px, drP.white, fixCenters(:, idx, b)', 2);

        % (g) Atualiza a tela para exibir a cruz de fixação
                        if mode == 1
                            Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
                            Screen('DrawTexture', dpP.window, bg); 
                            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                            lastPos = [-1 -1];
                        end
                        % disp('Vai calcular FSOnset')
                        FSonset = Screen('Flip', dpP.window, GetSecs() + .5*dpP.ifi);
                        FPonset = FSonset;
                    
        % iv. Inicia o registro da sessão
                        if debug == 0 && mode >= 2
                            Eyelink('StartRecording');
                            wasRecording = true;
                            Eyelink('Message',sprintf(prm.msg.on.trl{1}, i, tkP.nTrials));
                            Eyelink('Message',prm.msg.on.P1);
                        end
                    
        % v. Não avança de tela até que o olho (ou o cursor) esteja na
        %    cruz de fixação
                        fixCenter = fixCenters(:,idx, b);
                         if debug > 0 && mode >= 2
                            while true
                                [keyIsDown, ~, keyCode] = KbCheck;
                                KbReleaseWait;
                                if keyIsDown
                                    if keyCode(escapeKey)
                                        [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, txP, dpP, drP, prm, 'trial', debug, mode, targetOri(b));
                                    else
                                        fixAcquired = true;
                                    end
                                    break;
                                end
                                WaitSecs(0.02);
                            end
                        else
                            while true
                                check = false;
                                % Se exceder o tempo limite sem alterar
                                % fixAcquired, o while ~fixAcquired será
                                % reiniciado
                                if (GetSecs - FSonset) > prm.maxCrossDur
                                    break;
                                end
                                if mode >= 2
                                    damn = Eyelink('CheckRecording');
                                    if(damn ~= 0), break; end
            
                                    if Eyelink('NewFloatSampleAvailable') > 0
                                        evt = Eyelink('NewestFloatSample');
                                        x_gaze = evt.gx(tkP.Eye);
                                        y_gaze = evt.gy(tkP.Eye);
                                        check = true;
                                    end
                                elseif mode == 1
                                    [x_gaze, y_gaze, ~] = GetMouse(dpP.window);
                                    check = true;
                                    if any([x_gaze, y_gaze] ~= lastPos)
                                        Screen('DrawTexture', dpP.window, bg);
                                        Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
                                        Screen('Flip', dpP.window, GetSecs() + .5*dpP.ifi);
                                        lastPos = [x_gaze, y_gaze];
                                    end
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
                                        [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, txP, dpP, drP, prm, 'trial', debug, mode, targetOri(b));
                                        break;
                                    end
                                end
                            end
                        end
                        if ~fixAcquired && keepGoingTrials && ~restartTrial
                            if debug == 0 && mode >= 2
                                Eyelink('Message',prm.msg.err.P1);
                                Eyelink('StopRecording');
                                Eyelink('SetOfflineMode');
                            else
                                disp('Erro: sem fixação no tempo necessário')
                            end
                                
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
                                        if debug == 0 && mode >= 2
                                            Eyelink('Message',prm.msg.pse{5});
                                            EyelinkDoTrackerSetup(tkP.el);
                                            restartTrial = true;
                                        else
                                            disp('Recalibragem solicitada')
                                        end
                                        break;
                                    elseif keyCode(spaceKey)
                                        if debug == 0 && mode >= 2
%                                             Eyelink('Message', 'PRETRIAL_NO_RECALIBRATION');
                                        else
                                            disp('Prossegue sem recalibragem')
                                        end
                                        break;
                                    elseif keyCode(escapeKey)
                                        [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, txP, dpP, drP, prm, 'trial', debug, mode, targetOri(b));
                                        break;
                                    end
                                end
                                WaitSecs(0.02);
                            end
                        end
                    end
                    
    %% 6) Início Fase 2: tela de estímulos
                    Screen('Flip', dpP.window, GetSecs() + .5*dpP.ifi);
                    if keepGoingTrials && ~restartTrial
                        % Novamente, no modo cursor é usada textura de tela cheia
                        if mode == 1
                            Screen('Close', bg); clear auxWin bg;
                            bg = Screen('OpenOffscreenWindow', dpP.window, drP.grey, [], 64);
                            auxWin = bg;
                        else
                            auxWin = dpP.window;
                        end

                        alphas = ones(nStims, 1);

                        foragingDrawMain(auxWin, gaborTex, noiseTex, srcRects, dstRects, orientation(:, idx, b), txP, [repmat(alphas', [3,1]); ones(1, nStims)]);
            
            % (k) Atualiza a tela para exibir os estímulos
                        if mode == 1
                            Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
                            Screen('DrawTexture',   dpP.window, bg); 
                            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                            lastPos = [-1 -1];
                        end
                        trialOnset = Screen('Flip', dpP.window, GetSecs() + 0.5 * dpP.ifi);
            
            % vi. Registra o momento de início do trial
                        if debug == 0  && mode >= 2 && keepGoingTrials
                            Eyelink('Message',[sprintf('I %d ', uint64(trialOnset * 1000)) prm.msg.on.P2]);
                        end
                    else
                        trialOnset = GetSecs;
                    end
        
        % vii. O sujeito deve visitar exatamente Ns(i) estímulos
        %      antes de ocorrerem as modificações pré-sacádicas
                    counter = 0;                 % Número de estímulos visitados
                    flag = zeros(1, nStims);     % Quantas vezes cada estímulo foi visitado
                    stimTimes = nan(nStims,1);
                    currStim = 0;
                    fixStartTime = NaN;
                    if keepGoingTrials && ~restartTrial
                        runTrial = true;
                        while runTrial

                            WaitSecs(0.001);
                            tNow = GetSecs;
                            alphasPrev = alphas;
                            alphas = getGaborAlpha(tNow, stimTimes, prm);
                            alphaChanged = ~isequal(alphasPrev, alphas);
                            if alphaChanged
                                if mode == 1
                                    Screen('FillRect', auxWin, drP.grey);
                                end
                                foragingDrawMain(auxWin, gaborTex, noiseTex, srcRects, dstRects, orientation(:, idx, b), txP, [repmat(alphas', [3,1]); ones(1, nStims)])
                            end
                        
                            check = false;
                            if mode >= 2
                                damn = Eyelink('CheckRecording');
                                if(damn ~= 0), break; end
                        
                                if Eyelink('NewFloatSampleAvailable') > 0
                                    evt = Eyelink('NewestFloatSample');
                                    x_gaze = evt.gx(tkP.Eye);
                                    y_gaze = evt.gy(tkP.Eye);
                                    check = true;
                                end
                                if alphaChanged
                                    Screen('Flip', dpP.window, tNow + .5*dpP.ifi);
                                end
                            elseif mode == 1
                                [x_gaze, y_gaze, ~] = GetMouse(dpP.window);
                                check = true;
                                if any([x_gaze, y_gaze] ~= lastPos) || alphaChanged
                                    Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
                                    Screen('DrawTexture',   dpP.window, bg); 
                                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                                    Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
                                    Screen('Flip', dpP.window, tNow + .5*dpP.ifi);
                                    lastPos = [x_gaze, y_gaze];
                                end
                            end
                        
                            if check
                                % Se o olho (ou cursor) estiver próximo de um
                                % alvo não antes visto por tempo suficiente, 
                                % incrementa o contador
                                isCurrStim = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, idx, b)) <= minFixDist1;
                                currIdx = find(isCurrStim, 1);
                        
                                % Quando não está fixando num estímulo mas
                                % estava, a fixação é registrada se longa
                                if isempty(currIdx)
                                    if currStim ~= 0
                                        fixDur = tNow - fixStartTime;
                            
                                        if fixDur >= prm.minFixTime2
                                            if debug == 0 && mode >= 2
                                                Eyelink('Message',prm.msg.off.stm{1});
                                                disp('Fim de fixação boa em stim');
                                            end
                                            % Apenas salva na fila a primeira
                                            % fixação em um estímulo
                                            if flag(currStim) == 0
                                                counter = counter + 1;
                                                fprintf('Terminou a visita ao %d-ésimo estímulo\n', counter)
                                                auxFixQueue(counter) = fixDur;
                                                [P3On, tkP] = P3Onset5(tkP, prm, fixDur);
                                            end
                                            seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currStim; fixDur; 2]]; % Se quisesse registrar o comprimento de todas as fixações
                                            flag(currStim) = flag(currStim) + 1;
                                            % Registra como ruim a fixação se tiver sido muito curta
                                        else
                                            if debug == 0 && mode >= 2
                                                Eyelink('Message',prm.msg.off.stm{2}); 
                                                disp('Fim de fixação ruim em stim')
                                                WaitSecs(0.002);
                                            end
                                        end
                                    end
                                    % Como deixou de fixar, reseta as variáveis
                                    % de estado e início de fixação
                                    currStim = 0;
                                    % fprintf('currStim = %d\n', currStim);
                                    fixStartTime = NaN;
                                
                                % Entra no else quando há estímulo fixado
                                % (i.e., começa ou continua fixando)
                                else
                                    % Se logo antes não havia estímulo fixado,
                                    % ou havia um estímulo diferente
                                    % (pouco provável, requer estímulos 
                                    % próximos E amostragem baixa), temos
                                    % o início da fixação
                                    if currStim == 0 || currStim ~= currIdx

                                        % Adiciona mensagem de fim de fixação
                                        % para esse caso improvável
                                        if currStim ~= 0
                                            disp('Atenção: pulou de um estímulo a outro')
                                            fixDur = tNow - fixStartTime;
                                            if fixDur >= prm.minFixTime2
                                                if debug == 0 && mode >= 2
                                                    Eyelink('Message',prm.msg.off.stm{1}); 
                                                    disp('Fim de fixação boa em stim')
                                                end
                                                if flag(currStim) == 0
                                                    counter = counter + 1;
                                                    fprintf('Terminou a visita ao %d-ésimo estímulo\n', counter)
                                                    auxFixQueue(counter) = fixDur;
                                                    [P3On, tkP] = P3Onset5(tkP, prm, fixDur);
                                                end
                                                seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currStim; fixDur; 2]]; % Se quisesse registrar o comprimento de todas as fixações
                                                flag(currStim) = flag(currStim) + 1;
                                            else
                                                if debug == 0 && mode >= 2
                                                    Eyelink('Message',prm.msg.off.stm{2});                                                     
                                                    disp('Fim de fixação ruim em stim')
                                                end
                                            end
                                            WaitSecs(0.002);
                                        end

                                        currStim = currIdx;
                                        fprintf('currStim = %d\n', currStim);
                        
                                        if debug == 0 && mode >= 2
                                            if flag(currStim) == 0
                                                Eyelink('Message',prm.msg.on.stm{1});
                                            else
                                                Eyelink('Message',prm.msg.on.stm{2});
                                            end
                                        end
                                        fixStartTime = tNow;
                                        if flag(currStim) == 0
                                            stimTimes(currStim) = fixStartTime;
                                        end
                        
                                        %% IMPORTANTE: Se for iniciada a modTimes(b,i)-ésima fixação
                                        % diferente num estímulo, começa a contar o tempo de
                                        % atualização
                                        if flag(currStim) == 0 && counter == modTimes(b, idx) - 1
                                            preUpdateDeadline = fixStartTime + P3On;
                                            counter = counter+1;
                                            fprintf('Visita ao último %d-ésimo estímulo\n', counter)
                                             fprintf('Tempo pré-ruído rosa disponível: %.4f\n', preUpdateDeadline - GetSecs)
                                        end
                        
                                        % Se ainda estiver fixando o mesmo estímulo
                                        % (i.e., currStim == currIdx), não faz nada
                                    end
                                end
                            else
                            end
                        % () Se o tempo máximo tiver sido excedido, exibe uma tela especial
                            tAux = tNow - trialOnset;
                            if tAux > maxTrialDur
                                if debug == 0 && mode >= 2
                                    Eyelink('Message',prm.msg.err.P2);
                                    Eyelink('StopRecording');
                                    Eyelink('SetOfflineMode');
                                else
                                    disp('Erro: tempo de busca excedido')
                                end
                                restartTrial = true;
                                Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                                Screen('FrameOval', dpP.window, drP.red, dstRects, prm.pW2);
                                
                                Screen('Flip', dpP.window, tNow + .5*dpP.ifi);
                                WaitSecs(prm.fadeInDelay1);
                                
                                tStart = GetSecs;
                                while true
                                    tNow = GetSecs;
                                    alpha = min((tNow - tStart) / prm.fadeInDur1, 1);
                        
                                    Screen('FrameOval', dpP.window, drP.red, dstRects, prm.pW2);
                        
                                    Screen('FillRect', dpP.window, [overlayColor alpha*.75*drP.white], overlayRect);
                                    Screen('TextSize', dpP.window, prm.textSizeEnormous); Screen('TextStyle', dpP.window, 1);
                                    DrawFormattedText(dpP.window, 'TEMPO ESGOTADO', 'center', dpP.winRect(4)/2 - 25, [textColor1 alpha*drP.white]);
                                    if i < tkP.nTrials && retryCount(trialQueue(i)) < prm.maxRetries
                                        Screen('TextSize', dpP.window, prm.textSizeBigger); Screen('TextStyle', dpP.window, 0);
                                        DrawFormattedText(dpP.window, 'Por favor, tente novamente', 'center', dpP.winRect(4)/2 + 50, [textColor1 alpha*drP.white]);
                                    end
                                    Screen('TextSize', dpP.window, prm.textSizeBig);
                                    Screen('Flip', dpP.window);
                                
                                    if alpha >= 1, break; end
                                end
                                WaitSecs(prm.fadeInDelay1); tStart = GetSecs;
                                while true
                                    tNow = GetSecs;
                                    alpha = min((tNow - tStart) / prm.fadeInDur2, 1);
                        
                                    Screen('FrameOval', dpP.window, drP.red, dstRects, prm.pW2);
                        
                                    Screen('FillRect', dpP.window, [overlayColor .75*drP.white], overlayRect);
                                    Screen('TextSize', dpP.window, prm.textSizeEnormous); Screen('TextStyle', dpP.window, 1);
                                    DrawFormattedText(dpP.window, 'TEMPO ESGOTADO', 'center', dpP.winRect(4)/2 - 25, [textColor1 drP.white]);
                                    if i < tkP.nTrials && retryCount(trialQueue(i)) < prm.maxRetries
                                        Screen('TextSize', dpP.window, prm.textSizeBigger); Screen('TextStyle', dpP.window, 0);
                                        DrawFormattedText(dpP.window, 'Por favor, tente novamente', 'center', dpP.winRect(4)/2 + 50, [textColor1 drP.white]);
                                    end
                                    Screen('TextSize', dpP.window, prm.textSizeBig);
                                    DrawFormattedText(dpP.window, 'Pressione ESPAÇO para prosseguir', 'center', dpP.winRect(4)/2 + 110, [textColor2 alpha*drP.white]);
                                    Screen('Flip', dpP.window);
                                
                                    if alpha >= 1, break; end
                                end
                                Screen('TextSize', dpP.window, prm.textSizeNormal);
                        
                                while true
                                    [keyIsDown, ~, keyCode] = KbCheck;
                                    KbReleaseWait;
                                    if keyIsDown
                                        if keyCode(spaceKey)
                                            break;
                                        elseif keyCode(escapeKey)
                                            [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, txP, dpP, drP, prm, 'trial', debug, mode, targetOri(b));
                                            break;
                                        end
                                    end
                                end
                            end
                            runTrial = (counter < modTimes(b,idx)) && ~restartTrial;
                        end
                    end
                    fprintf('currStim final = %d\n', currStim);
        
    
    %% 7) Início Fase 3: tela com ruído rosa
        % (l) Desenha ruído orientado em todos os estímulos menos o 
        %     atual -- as linhas comentadas servem para mudar apenas os 
        %     não vistos. Atrasa a apresentação para o estímulo durar
        %     medFixTime segundo antes de o ruído rosa substituí-lo
                    if mode == 1
                        auxWin  = bg;
                    else
                        auxWin = dpP.window;
                    end
                    if ~restartTrial && keepGoingTrials
                        blinkIdx = setdiff(1:(nStims), currIdx);
                        fprintf('Tempo permitido de fixação antes do rosa: %.4f\n', P3On)

                        
                        P3Frames = round(P3On / dpP.ifi);
                        TargetP3Onset = fixStartTime + (P3Frames - 0.5) * dpP.ifi;
                        
                        % O update ocorre 1 frame antes do início da fase 3
                        preUpdateDeadline = TargetP3Onset - dpP.ifi; 
                        tNow = GetSecs;
                        while tNow < preUpdateDeadline
                            WaitSecs(0.001);
                            tNow = GetSecs;
                            alphasPrev = alphas;
                            alphas = getGaborAlpha(tNow, stimTimes, prm);
                            alphaChanged = ~isequal(alphasPrev, alphas);
                            if alphaChanged && check
                                if mode == 1
                                    Screen('FillRect', auxWin, drP.grey);
                                end
                                foragingDrawMain(auxWin, gaborTex, noiseTex, srcRects, dstRects, orientation(:, idx, b), txP, [repmat(alphas', [3,1]); ones(1, nStims)])
                            end

                            check = false;
                            if mode >= 2
                                damn = Eyelink('CheckRecording');
                                if(damn ~= 0), break; end
                        
                                if Eyelink('NewFloatSampleAvailable') > 0
                                    evt = Eyelink('NewestFloatSample');
                                    x_gaze = evt.gx(tkP.Eye);
                                    y_gaze = evt.gy(tkP.Eye);
                                    check = true;
                                end
                                if alphaChanged && check
                                    Screen('Flip', dpP.window, tNow + .5*dpP.ifi);
                                end
                            elseif mode == 1
                                [x_gaze, y_gaze, ~] = GetMouse(dpP.window);
                                check = true;
                                if any([x_gaze, y_gaze] ~= lastPos) || alphaChanged
                                    Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
                                    Screen('DrawTexture',   dpP.window, bg); 
                                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                                    Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
                                    Screen('Flip', dpP.window, tNow + .5*dpP.ifi);
                                    lastPos = [x_gaze, y_gaze];
                                end
                            end
                        end
                        if mode >= 3
                            Screen('BlendFunction', auxWin, GL_ONE, GL_ONE);
                            Screen('DrawTextures', auxWin, oriPinkTex, srcRects(:,blinkIdx), dstRects(:,blinkIdx), orientation(blinkIdx, idx, b), [], [], [], []);
    
                            if ~isempty(currIdx)
                                Screen('BlendFunction', auxWin, GL_ONE, GL_ONE);
                                Screen('DrawTextures', auxWin, gaborTex, [], dstRects(:,currIdx), orientation(currIdx, idx, b));
                                Screen('DrawTextures', auxWin, noiseTex, srcRects(:,currIdx), dstRects(:,currIdx), orientation(currIdx, idx, b));
                                
                                Screen('BlendFunction', auxWin, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                                Screen('DrawTextures', auxWin, txP.blob.tex, [], dstRects(:, currIdx), orientation(currIdx, idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
                            else
                                Screen('BlendFunction', auxWin, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                            end
                        
                % (j) Desenha a abertura gaussiana
                            Screen('DrawTextures', auxWin, txP.blob.tex, [], dstRects(:, blinkIdx), orientation(blinkIdx, idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
                            Screen('Close', oriPinkTex); Screen('Close', gaborTex);
                            Screen('BlendFunction', auxWin, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        end
                        if mode >= 3
                            updateStimOnset = Screen('Flip', dpP.window, TargetP3Onset);
                        else
                            WaitSecs('UntilTime', TargetP3Onset);
                            updateStimOnset = GetSecs;
                        end

                        if debug == 0 && mode >= 2, Eyelink('Message',[sprintf('I %d ', uint64(updateStimOnset * 1000)) prm.msg.on.P3]); end
                        
                        % A fase 3 é encerrada se o estímulo fica tempo
                        % demais na tela ou quando o olho sai do último
                        % estímulo
                        fprintf('currStim último fixado (agora em P3) = %d\n', currStim);

                        P3Deadline = updateStimOnset + (prm.pinkNoiseNFrames - 1)*dpP.ifi;
                        earlyStop = false;
                        while true
                            tNow = GetSecs;
                            if tNow > P3Deadline
                                P3Dur = tNow - updateStimOnset;
                                fixDur = tNow - fixStartTime;
                                fprintf('Fim ruído rosa por duração: %.4f\n', P3Dur)
                                if debug == 0 && mode >= 2
                                    Eyelink('Message',prm.msg.off.stm{3});
                                    disp('Fim de fixação P3 em stim')
                                end
                                break;
                            end
                            check = false;
                            if mode >= 2
                                damn = Eyelink('CheckRecording');
                                if(damn ~= 0), break; end
                
                                if Eyelink('NewFloatSampleAvailable') > 0
                                    evt = Eyelink('NewestFloatSample');
                                    x_gaze = evt.gx(tkP.Eye);
                                    y_gaze = evt.gy(tkP.Eye);
                                    check = true;
                                end
                            elseif mode == 1
                                [x_gaze, y_gaze, ~] = GetMouse(dpP.window);
                                check = true;
                                if any([x_gaze, y_gaze] ~= lastPos)
                                    Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
                                    Screen('DrawTexture', dpP.window, bg);
                                    Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
                                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                                    Screen('Flip', dpP.window, tNow + .5*dpP.ifi);
                                    lastPos = [x_gaze, y_gaze];
                                end
                            end
                            if check
                                isCurrStim = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, idx, b)) <= minFixDist1;
                                if isCurrStim(currStim) == 0
                                    P3Dur = tNow - updateStimOnset;
                                    fixDur = tNow - fixStartTime;
                                    fprintf('Fim ruído rosa por dispersão: %.4f\n', P3Dur)
                                    earlyStop = true;
                                    if debug == 0 && mode >= 2
                                        Eyelink('Message',prm.msg.off.stm{1});
                                        disp('Fim de fixação P3 em stim')
                                    end
                                    break;
                                end
                            end
                            WaitSecs(.0005);
                        end
                        isP3earlyStop(b, idx) = earlyStop;
                        % Como saio do loop da fase 2 assim que inicio a
                        % última fixação, tenho que adicionar neste momento
                        % a fixação iniciada lá
                        seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currStim; fixDur; 3]];
                        flag(currStim) = flag(currStim) + 1;

                        foragingDrawPedestal(dpP.window, noiseTex, srcRects, dstRects, orientation(:, idx, b), txP);
                        Screen('Close', noiseTex);
%                         Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                        if earlyStop || mode == 1
                            P3Flip = GetSecs + 0.5 * dpP.ifi;
                        else
                            P3Flip = updateStimOnset + (prm.pinkNoiseNFrames - 0.5) * dpP.ifi;
                        end
                        updateStimOffset = Screen('Flip', dpP.window, P3Flip);

                        if debug == 0 && mode >= 2,  Eyelink('Message',[sprintf('I %d ', uint64(updateStimOffset * 1000)) prm.msg.off.P3]); end

                        P3Dur = updateStimOffset - updateStimOnset;
                        fprintf('Tempo total de ruído rosa: %.4f\n', P3Dur)

                        seenIdx = find(flag ~= 0);
                        notSeenIdx = find(flag == 0);
                        
            
            %% ix. Verifica se há alguma fixação de duração mínima em estímulo 
            %     numa janela pós-modificação
                        fixOnset = updateStimOffset; currIdx = [];
                        initialIdx = 0; % Para rastrear o primeiro ROI fixado

                        if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.on.PM); end
            
                        if debug == 0 || mode == 1
                            maxDurReached = true;
                            tNow = GetSecs;
                            while tNow - updateStimOffset < prm.postModDur
                                check = false;
                                if mode >= 2
                                        
                                    damn = Eyelink('CheckRecording');
                                    if(damn ~= 0), break; end
                        
                                    if Eyelink('NewFloatSampleAvailable') > 0
                                        evt = Eyelink('NewestFloatSample');
                                        x_gaze = evt.gx(tkP.Eye);
                                        y_gaze = evt.gy(tkP.Eye);
                                        check = true;
                                    end
                                elseif mode == 1
                                    [x_gaze, y_gaze, ~] = GetMouse(dpP.window);
                                    check = true;
                                end
                        
                                if check
                                    isCurrStim = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, idx, b)) <= minFixDist3;
                        
                                    % Se currIdx estava com valor inicial e isCurrStim não é  
                                    % totalmente nulo, começou uma fixação 
                                    if isempty(currIdx)
                                        if any(isCurrStim)
                                            currIdx = find(isCurrStim, 1);
                                            % Guarda a primeira fixação
                                            if initialIdx == 0
                                                initialIdx = currIdx; 
                                            end
                        
                                            fixOnset = tNow;
                                            if debug == 0 && mode >= 2
                                                if flag(currIdx) == 0
                                                    Eyelink('Message',prm.msg.on.stm{1});
                                                else
                                                    Eyelink('Message',prm.msg.on.stm{2});
                                                end
                                            end
                                        end
                                    % Uma vez que currIdx é não nulo, fica verificando se a
                                    % fixação saiu dele, o que acontece se isCurrStim(currIdx)
                                    % voltar a ser nulo. Considera como fixação apenas se for
                                    % suficientemente longa
                                    else
                                        if ~isCurrStim(currIdx)
                                            fixDur = tNow - fixOnset;
                                            % Se a fixação não havia escapado do alvo ainda, é só
                                            % uma questão de latência e não há porque encerrar a PM
                                            if currIdx == initialIdx
                                                disp(['Saiu do currIdx ' num2str(currIdx) ', temos tempo sobrando em PM']);
                                                if debug == 0 && mode >= 2
                                                    Eyelink('Message',prm.msg.off.stm{1}); 
                                                    disp('Fim de fixação PM em stim')
                                                end
                                            else
                                                if fixDur >= prm.minFixTime3
                                                    seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currIdx; fixDur; 4]];
                                                    disp(['Trial ' num2str(idx) ': Visitou o alvo ' num2str(currIdx) ' pós-modificação']);
                                                    if debug == 0 && mode >= 2
                                                        Eyelink('Message',prm.msg.off.stm{1}); 
                                                        disp('Fim de fixação PM em stim')
                                                    end
                                                    maxDurReached = false;
                                                    break
                                                else
                                                    if debug == 0 && mode >= 2
                                                        Eyelink('Message',prm.msg.off.stm{2});
                                                        disp('Fim de fixação ruim PM em stim')
                                                    end
                                                end
                                            end
                                            % Se a fixação não foi longa o suficiente, pode ser
                                            % que estava apenas passando pelo estímulo para 
                                            % fovear outro, então recomeça a busca
                                            currIdx = [];
                                        else
                                            fixDur = tNow - fixOnset;
                                            % Se ainda está na fixação inicial e excedeu temp de tolerância, esquece 
                                            if (currIdx == initialIdx) && (fixDur >= prm.maxTolPM)
                                                seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currIdx; fixDur; 4]];
                                                disp([' Sem movimento ocular em PM. Permaneceu no alvo ' num2str(currIdx)]);
                                                if debug == 0 && mode >= 2
                                                    Eyelink('Message', prm.msg.off.stm{1}); 
                                                    disp('Fim de fixação PM por inércia (maxTolPM)')
                                                end
                                                maxDurReached = false;
                                                break;
                                            end
                                        end
                                    end
                        
                        
                                    if mode == 1
                                        [x_gaze, y_gaze, ~] = GetMouse(dpP.window);
                                        if any([x_gaze, y_gaze] ~= lastPos)
                        %                                             if tNow - updateStimOffset <= prm.blobPMDur
                                            Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                        %                                             end
                                            Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
                                            Screen('Flip', dpP.window, tNow + .5*dpP.ifi);
                                            lastPos = [x_gaze, y_gaze];
                                        end
                                    end
                                end
                                WaitSecs(.0005);
                                tNow = GetSecs;
                            end
                            if maxDurReached && ~isempty(currIdx)
                                seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currIdx; prm.postModDur; 4]];
                                disp(['Trial ' num2str(idx) ': Visitou o alvo ' num2str(currIdx) ' pós-modificação (fim forçado)']);
                                if debug == 0 && mode >= 2
                                    Eyelink('Message',prm.msg.off.stm{4}); 
                                    disp('Fim de fixação PM em stim')
                                end
                            elseif ~isempty(currIdx)
                                 WaitSecs(prm.postModDur - (tNow - updateStimOffset));
                            end
                        end
                        fprintf('seenStimsQueue final: '); disp(seenStimsQueue{b,i}(1,:)); fprintf('\n')

                        if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.PM); end
%                         Screen('Flip', dpP.window);
                
                            if ~isempty(currIdx)
                                notSeenIdx(notSeenIdx == currIdx) = [];
    
                                % Se não tiver mexido os olhos, não será feita
                                % nenhuma pergunta sobre o estímulo fixado
                                if currIdx == currStim
                                    seenIdx(seenIdx == currIdx) = [];
                                    currIdx = [];
                                    disp('Não chegou noutro estímulo durante ruído rosa')
                                end
                            end
                        
    %% 8) Início Fase 4: reportar em quais posições havia alvos
            % x. Desenha placeholders para os estímulos aleatoriamente, obedecendo
            %    a identificação deles como visitados, atual e não visitados
                        if mode == 1
                            Screen('Close', bg); clear auxWin bg;
                        end

                        % % Se houver fixado em algum estímulo em PM...
                        skipP4 = true;
                        if ~isempty(currIdx)
                            skipP4 = false;
                            fprintf('Há currIdx... ')
                            fprintf('então será perguntado sobre\n');
                        else
%                             nPre = nStimsToReport(1, idx, b); nPost = nStimsToReport(3, idx, b);
                            fprintf('Não há currIdx... ')
                            if mode <= 2
                                if rand < prm.seenNotSeenRatio/(1+prm.seenNotSeenRatio)
                                    fprintf('... então pré-s. vira visto\n')
                                    nPre = nStimsToReport(1, idx, b) + nStimsToReport(2, idx, b);
                                    nPost = nStimsToReport(3, idx, b);
                                else
                                    fprintf('... então pré-s. vira não visto\n')
                                    nPre = nStimsToReport(1, idx, b);
                                    nPost = nStimsToReport(3, idx, b) + nStimsToReport(2, idx, b);
                                end
                                nStimsToReport(1, idx, b) = nPre;
                                nStimsToReport(2, idx, b) = 0;
                                nStimsToReport(3, idx, b) = nPost;
                            
                                if numel(seenIdx) < nPre
                                    dif = nPre - numel(seenIdx);
                                    nStimsToReport(1, idx, b) = numel(seenIdx);
                                    nStimsToReport(3, idx, b) = nPost + dif;
                                end
                                if numel(notSeenIdx) < nPost
                                    dif = nPost - numel(notSeenIdx);
                                    nStimsToReport(3, idx, b) = numel(notSeenIdx);
                                    nStimsToReport(1, idx, b) = nPre + dif;
                                end
                            end
                        end
%                         disp(nStimsToReport(:, idx, b))
                        if sum(nStimsToReport(:, idx, b)) ~= 3
                            disp('ERRO')
                        end

                        % Para o treino de forrageamento, não faz sentido
                        % punir com o vermelhinho, até porque ainda está 
                        % sendo feito o ajuste dos parâmetros temporais
                        if mode <= 2, skipP4 = false; end

                        if skipP4
                            warningFlip(dpP.window, stimCenters(:, :, idx, b), resizeRect(dstRects, .25), currStim, txP.gabor.size_px, drP.allPW, drP.darkRed);
                        else
                            isTargetAnswer = nan(1,nStims); allColors2 = drP.allColors;
                            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
                            % Os demais pontos que não o fixado antes da atual
                            nbhd = NeighborsOrder(stimCenters(:, :, idx, b), currStim, currIdx);
    
                            % O do passado não importa de onde eu pergunto,
                            % desde que não seja a última fixação nem o pós
                            % sacádico
                            isSaccSeen(b, idx) = any(ismember(currIdx, seenIdx));
                            auxIdx = setdiff(seenIdx, [currStim currIdx]);
                            seenAux = datasample(auxIdx, min(length(auxIdx), nStimsToReport(1, idx, b)), 'Replace', false);

                            currAux = []; 
                            if nStimsToReport(2, idx, b) == 1, currAux = currIdx; end

                            % O não visto eu pergunto o mais próximo de
                            % currStim que não seja currIdx e nem já visto. Mas
                            % se a pessoa nao mover o olho, como puniçao,
                            % perguntamos sobre os mais distantes
                            if isempty(currIdx)
                                nbhd = flip(nbhd);
                            end
    
                            notSeenNbhd = setdiff(nbhd, currIdx, 'stable');
                            notSeenNbhd = setdiff(notSeenNbhd, seenIdx, 'stable');
                            notSeenAux = notSeenNbhd(1:min(numel(notSeenNbhd),nStimsToReport(3, idx, b)));
    
                            % notSeenAux = datasample(auxNotSeenIdx, min(length(auxNotSeenIdx), nStimsToReport(3, idx, b)), 'Replace', false);
                            % Tanto seenAux como notSeenAux devem ser linhas
                            % para concatenar
                            fprintf('Vistos: '); fprintf(num2str(seenAux));
                            fprintf('\nPré-s: '); fprintf(num2str(currAux));
                            fprintf('\nNão vistos: '); fprintf(num2str(notSeenAux));
                            fprintf('\n');
                            orderToReportStimsCell = {seenAux, currAux, notSeenAux};
                            orderToReportStims = [orderToReportStimsCell{orderToReportSets(1, idx, b)} orderToReportStimsCell{orderToReportSets(2, idx, b)} orderToReportStimsCell{orderToReportSets(3, idx, b)}];
                            if numel(orderToReportStims) < 2 || numel(orderToReportStims) > 3
                                fprintf('Algo de errado: tem que reportar: %d', numel(orderToReportStims));
                            end
                            orderRemapped = [];
                            for auxIdx = 1:3
                                orderRemapped = [orderRemapped orderToReportMap(orderToReportSets(auxIdx, idx, b))*ones(1, length(orderToReportStimsCell{orderToReportSets(auxIdx, idx, b)}))]; %#ok<AGROW> 
                            end
    
                             if debug == 0  && mode >= 2
                                Eyelink('Message',prm.msg.on.P4);
                             end
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
                            
                                    isTargetAnswer(orderToReportStims(j)) = currTarget;
                                    Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                                    foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, rectColors, isTargetAnswer, targetOri(b), rectPW);
                            
                                end
                            end
                            
                            if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.P4); end
    
                            feedback = (orientation(:,idx,b) == targetOri(b))' + isTargetAnswer;
                            feedback(rem(feedback,2) == 0) = 2; feedback(rem(feedback,2) == 1) = 0;
                            feedback = feedback/2;
    
                            trialFeedback{b, i} = [orderToReportStims; orderRemapped; feedback(orderToReportStims)];
                            
                            if mode <= 3
                                Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                                foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, isTargetAnswer, targetOri(b), drP.allPW, feedback, drP.red, drP.green);
                            else
                                Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                                foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, isTargetAnswer, targetOri(b), drP.allPW);
                            end
                            
                            trialOrder(2, i, b) = 1;
    
                            i = i + 1;
                            trialIdxUp = true;
                        end

                        WaitSecs(.5);

                        restartTrial = restartTrial | skipP4;
                    % A condição é keepGoingTrials = false ou restartTrial = true
                    end

                    if keepGoingTrials == false || restartTrial == true
                        
                        if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.err.trl{1}); end
                        retryCount(trialQueue(i)) = retryCount(trialQueue(i)) + 1;
                        % Só atualizo a fila de trials se eles continuarem
                        % sendo úteis, i.e., se a única solicitação foi
                        % reiniciar o trial
                        if keepGoingTrials
                            if retryCount(trialQueue(i)) > prm.maxRetries
                                if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.err.trl{2}); end
                                warning('Trial %d excede o máximo de repetições. Prosseguindo', trialQueue(i));
                                i = i + 1;
                                trialIdxUp = true;
                            else
                                % Faz com que o trial não terminado sempre vá para
                                % o fim da fila
                                trialQueue(i:end) = [trialQueue(i+1:end) trialQueue(i)];
                            end
                        end
                     end
    
            % xi. Interrompe o registro, pois ou a tela será atualizada ou
            %     acabaram os trials
                    if debug == 0 && mode >= 2% && keepGoingTrials
                        trialOffset = GetSecs;
                        Eyelink('Message', sprintf(prm.msg.off.trl{1}, i - trialIdxUp, tkP.nTrials));
                        Screen('Flip', dpP.window, trialOffset + .5*dpP.ifi);
                        WaitSecs(0.1);
                        Eyelink('SetOfflineMode');
                        Eyelink('StopRecording');
                    else 
                        Screen('Flip', dpP.window, GetSecs + .5*dpP.ifi);
                        WaitSecs(0.5);
                    end
                end
%                 if keepGoingTrials
                % Veja que como restartBlock = true é sempre acompanhado de
                % keepGoingTrials = false, ao reiniciar um bloco vai haver
                % tanto erro do último trial como do bloco
                if restartBlock
                    if debug == 0 && mode >= 2
                        Eyelink('Message',prm.msg.err.blk);
                        Eyelink('Message',sprintf(prm.msg.off.blk{1}, b, tkP.nBlocks)); 
                    end
                    restartBlock = false;
                else
                    blockOffset = GetSecs; %#ok<NASGU>
                    if debug == 0 && mode >= 2, Eyelink('Message',sprintf(prm.msg.off.blk{1}, b, tkP.nBlocks)); end
                    b = b+1;
                end
%                 end
            end

            if b == tkP.nBlocks+1 && keepGoingBlocks
                if mode >= 2, Eyelink('Message',sprintf(prm.msg.off.ses{1}, suffix)); end
                blocksCompleted = true;
            else
                if mode >= 2, Eyelink('Message',sprintf(prm.msg.off.ses{2}, suffix)); end
            end

            if mode >= 2
                tkS(mode, 2) = blocksCompleted;
            end
            if debug == 0 && mode >= 2

                results.fixCenters = fixCenters;
                results.stimCenters = stimCenters;
                results.orientation = orientation;
                results.nTs = nTs;
                results.targetOri = targetOri;
                results.modTimes = modTimes;
                results.nStimsToReport = nStimsToReport;
                results.orderToReportSets = orderToReportSets;
                results.trialOrder = trialOrder;
                results.isP3earlyStop = isP3earlyStop;
                results.isSaccSeen    = isSaccSeen;
                results.trialFeedback = trialFeedback;
                results.seenStimsQueue = seenStimsQueue;
                results.nSnbhd = nSnbhd;
            else
                results = [];
            end
        catch
            if ~exist('fixCenters', 'var'),   fixCenters = []; end
            if ~exist('stimCenters', 'var'),  stimCenters = []; end
            if ~exist('orientation', 'var'),  orientation = []; end
            if ~exist('nTs', 'var'),       nTs = []; end
            if ~exist('targetOri', 'var'), targetOri = []; end
            if ~exist('modTimes', 'var'),  modTimes = []; end
            if ~exist('nStimsToReport', 'var'),     nStimsToReport = []; end
            if ~exist('orderToReportSets', 'var'),  orderToReportSets = []; end
            if ~exist('trialOrder', 'var'),      trialOrder = []; end
            if ~exist('isP3earlyStop', 'var'),   isP3earlyStop = []; end
            if ~exist('isSaccSeen', 'var'),      isSaccSeen = []; end
            if ~exist('trialFeedback', 'var'),   trialFeedback = []; end
            if ~exist('seenStimsQueue', 'var'),  seenStimsQueue = []; end

            results.fixCenters = fixCenters;
            results.stimCenters = stimCenters;
            results.orientation = orientation;
            results.nTs = nTs;
            results.targetOri = targetOri;
            results.modTimes = modTimes;
            results.nStimsToReport = nStimsToReport;
            results.orderToReportSets = orderToReportSets;
            results.trialOrder = trialOrder;
            results.isP3earlyStop = isP3earlyStop;
            results.isSaccSeen = isSaccSeen;
            results.trialFeedback = trialFeedback;
            results.seenStimsQueue = seenStimsQueue;
            
            if debug ~= 0, tkS(:) = 0; end
            tkP = foragingSave(tkS, 2, prm, dpP, drP, tkP, txP, results);
            psychrethrow(psychlasterror);
        end

    if mode <= 3
        tkP.nTrainBlocks = tkP.nBlocks;
        tkP.nBlocks = nBlocks;
        tkP.nTrainTrials = tkP.nTrials;
        tkP.nTrials = nTrials;
        if debug == 0, HideCursor(dpP.window); end

        % Se estiver no modo treino, ajusta a duração do ruído rosa
        % conforme as fixações salvas na fila
        if mode == 2 || mode == 3
            disp('Ajuste de pinkNoiseDur')
            shortFixLimit = prctile(tkP.fixQueue, prm.shortFixPerc);
            tkP.pinkNoiseDur = min(prm.pinkNoiseDur, shortFixLimit);
        end
    end
end
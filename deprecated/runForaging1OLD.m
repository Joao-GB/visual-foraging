function [tkP, tkS, results] = runForaging1(tkP, dpP, drP, txP, prm, debug, mode, tkS)
        % O modo default é experiment
        if nargin < 7, mode = 'experiment'; end
        tic;


        Screen('Flip', dpP.window);

        if strcmp(mode, 'experiment')
            fakeAux = startFake(tkP, dpP, drP, prm, 'start', [], []);
        end

        if isfield(tkP,'pinkNoiseDur'), prm.pinkNoiseDur = tkP.pinkNoiseDur; end

        modeMap = containers.Map({'cursor', 'training', 'experiment'}, 1:3);
        mode = modeMap(mode);

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
    
    %% 1) Se não estiver no modo experimental, ajusta o número de trials e blocos
        if mode < 3
            nBlocks = tkP.nBlocks;
            nTrials = tkP.nTrials;
            if mode == 1
                prm.minFixTime2 = prm.cursorMinFixTime2;
                prm.medFixTime2 = prm.cursorMedFixTime2;
                prm.minFixTime3 = prm.cursorMinFixTime3;
                prm.postModDur  = prm.cursorPostModDur;
                prm.pinkNoiseDur= prm.cursorPinkNoiseDur;
                tkP.nBlocks = 1;
                tkP.nTrials = prm.nTrialsTrain;
            % No treino, todas as orientações são alvo uma vez
            elseif mode == 2
                L = numel(prm.allOri);
                tkP.nBlocks = L;
                tkP.nTrials = max(prm.nTrialsTrain, ceil(1*prm.fixTimeQueueSize/L));
            end
        end
        fprintf('Tempo inicial: %.5f\n', toc)
        tic;
    %% 2) Define as distribuições das condições dos trials e dos blocos
        nTrialsBuffered = tkP.nTrials + prm.nBufferTrials;
        [nTs, nStims, targetOri, modTimes, nStimsToReport, orderToReportSets] = getForagingDistributions1(tkP.nStims, tkP.nMaxFix, nTrialsBuffered, tkP.nBlocks, prm);
        if mode == 2,  targetOri = prm.allOri(randperm(tkP.nBlocks)); end

        % Redefine os gabores conforme a quantidade de estímulos
        txP.gabor     = foragingGabor(dpP.window, dpP.screenRes, nStims, prm, dpP.monitorW_mm);
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
        % Versão anterior para escolher distratores com base nas demais
        % orientações
%             auxAllOri = prm.allOri;
%             auxAllOri(auxAllOri == targetOri(b)) = [];
%             orientation(:,:,b) = auxAllOri(randi(length(auxAllOri), tkP.nStims, nTrialsBuffered));
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
        minFixDist1 = (txP.gabor.size_px/2)*prm.fixDistFactor1;

        fprintf('Tempo rects a centers: %.5f\n', toc);
        tic;
        for b=1:tkP.nBlocks
            for i=1:nTrialsBuffered
        % (d) Matrizes com os centros dos estímulos (i.e., centros dos dstRects)
        %     e das cruzes de fixação
%                 [currFixCenter, currStimCenter] = getStimLocations1(dpP.ellipseProps, txP.gabor.size_px, minDist_px);
                [currFixCenter, currStimCenter, rMax] = getStimLocations2(dpP.ellipseProps, nStims, minDist_px, false);
                fixCenters(:, i, b)     = currFixCenter;
                stimCenters(:, :, i, b) = currStimCenter;
        % (e) Matriz com orientações tem os alvos adicionados
                orientation(randperm(nStims, nTs(b, i)), i, b) = targetOri(b);
            end
        end

        % (f) Distância mínima para considerar fixação em alvo
        % pós-modificação (raio, não diâmetro)
        minFixDist3 = (txP.gabor.size_px/2)*prm.fixDistFactor3;

        auxFixQueue = zeros(1, nStims);
        if strcmp(mode, 'experiment')
            endFake(fakeAux, dpP, drP, prm);
        end
        fprintf('Tempo preencher locais estímulos e orientações de alvos: %.5f\n', toc);
    %% 4) Início dos blocos e trials
        fprintf('----Início da sessão----\n')
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

            if mode >= 2
                Eyelink('Message',sprintf(prm.msg.on.ses{1}, suffix));
                EyelinkDoTrackerSetup(tkP.el);
            end

            while b <= tkP.nBlocks && keepGoingBlocks
                trialOrder(:,:,b) = 0;
                
                % Para as telas de erro
                overlayRect = CenterRectOnPointd([dpP.winRect(1:3) dpP.winRect(4)/2.5], dpP.winRect(3)/2, dpP.winRect(4)/2);
                overlayColor = repmat(drP.whiteGrey, 1, 3);
                textColor1   = repmat(drP.black, 1, 3);
                textColor2   = repmat(.5*(drP.grey+drP.black), 1, 3);

        % (a) Mensagem com orientação a ser procurada
                % É repetida se o sujeito pausar a tarefa e voltar para o
                % bloco (seja recalibrando ou não)
                repeatMessage = true;
                while repeatMessage && keepGoingBlocks
                    Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
                    Screen('DrawTexture', dpP.window, txP.exampleNoise.tex, [], [], targetOri(b), [], [], [], []);
                    Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
                    Screen('DrawTexture', dpP.window, txP.exampleGabor.tex, [], [], targetOri(b), [], [], [], [], [], txP.exampleGabor.props);
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
                
                if debug == 0 && mode >= 2 && keepGoingBlocks
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
                    if modTimes(b, idx) >= tkP.nStims
                        fprintf('Algo errado!')
                    end

        % (b) Obtém a mediana do vetor de tempos de fixação
                    % auxN = [];
                    % if i > 1
                    %     auxN = size(seenStimsQueue{b, i-1},2);
                    % end
                    % P3On = P3Onset1(tkP, prm, auxN);
                    med  = median(tkP.fixQueue);
                    medFixTime = median(tkP.fixQueue);
                    if mode == 1
                        maxTrialDur = prm.cursorMaxTrialDurFactor*tkP.nStims*medFixTime;
                    else
                        maxTrialDur = prm.maxTrialDurFactor*tkP.nStims*medFixTime;
                    end
        % (c) Cria os retângulos de destino com base nas coordenadas dos
        %     centros dos estímulos
                    dstRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', stimCenters(1,:, idx, b), stimCenters(2,:, idx, b));
                    
                
                    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    if debug == 0 && mode >= 2
        % i. Deixa de registrar até StartRecording
                        if wasRecording, Eyelink('StopRecording'); end
                        Eyelink('SetOfflineMode');
                        WaitSecs(.1);
        % ii. Faz drift correction
                        EyelinkDoDriftCorrection(tkP.el);
                        WaitSecs(.1);
    
        % iii. Desenha, na tela do Host PC, linhas indicando a orientação de cada estímulo
                        Eyelink('Command', 'clear_screen 0');
                        for j=1:size(dstRects, 2)
                            hostRect = dstRects(:,j);
                            cX = mean(hostRect([1 3]));
                            cY = mean(hostRect([2 4]));
                            hostRect = round(hostRect);
                            radius = (hostRect(3)-hostRect(1))/2;
                            theta = orientation(j, idx, b); lineLen = radius * 0.8;
                            dx = lineLen * sind(theta); dy = lineLen * cosd(theta);
                            Eyelink('Command', 'draw_line %d %d %d %d 12', round(cX-dx), round(cY-dy), round(cX+dx), round(cY+dy));
                        end
                    end
                    
    %% 5) Início Fase 1: tela de fixação
                    % Se estiver no modo de cursor, cria uma textura em janela
                    % offscreen sobre a qual será desenhado o símbolo do cursor
                    if mode == 1
                        bg = Screen('OpenOffscreenWindow', dpP.window, drP.grey);
                        auxWin  = bg;
                    else
                        auxWin = dpP.window;
                    end

        % (d) Cria as texturas de ruído (pois não faz sentido armazená-las), sem desenhá-las
                    auxNoiseMatrix = butterFilter(pinkNoise(txP.gabor.size_px, (nStims)*txP.gabor.size_px), txP.noiseLoCutFreq, txP.noiseHiCutFreq);
                    noiseMatrix(:,:) = prm.noiseAlpha*(prm.noiseContrast*rescale(auxNoiseMatrix, -prm.noiseAmplitude, prm.noiseAmplitude))+drP.grey;
                    noiseTex = Screen('MakeTexture', dpP.window, noiseMatrix, [], [], [], [], []);
        
        % (e) Cria as texturas dos ruídos com orientação, sem desenhá-las
                    for j=1:(nStims)
                        colRange = ((j-1)*txP.gabor.size_px+1):(j*txP.gabor.size_px);
                        oriPinkMatrix(:,colRange) = ApplyOriFilter(txP.oriFilter', txP.OFsize, auxNoiseMatrix(:,colRange));
                    end
                    oriPinkTex = Screen('MakeTexture', dpP.window, oriPinkMatrix);
    
        % (f) Desenha a cruz de fixação
                    xFix = [-crossSize_px/2 crossSize_px/2 0 0];
                    yFix = [0 0 -crossSize_px/2 crossSize_px/2];
                    fixCoords = [xFix; yFix];
                    fixAcquired = false;
                    while ~fixAcquired && keepGoingTrials && ~restartTrial
                        Screen('BlendFunction', auxWin, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        Screen('DrawLines', auxWin, fixCoords, prm.lineWidth_px, drP.white, fixCenters(:, idx, b)', 2);
        
        % (g) Atualiza a tela para exibir a cruz de fixação
                        if mode == 1, Screen('DrawTexture', dpP.window, bg); lastPos = [-1 -1]; end
                        % disp('Vai calcular FSOnset')
                        FSonset = Screen('Flip', dpP.window);
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
                         if debug > 0 && mode > 1
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
                                if mode > 1
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
                                        Screen('Flip', dpP.window);
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
                    vbl = Screen('Flip', dpP.window);
                    if keepGoingTrials && ~restartTrial
                        % Novamente, no modo cursor é usada textura de tela cheia
                        if mode == 1
                            Screen('Close', bg); clear auxWin bg;
                            bg = Screen('OpenOffscreenWindow', dpP.window, drP.grey);
                            auxWin  = bg;
                        else
                            auxWin = dpP.window;
                        end
            
            % (h) Desenha os ruídos, puramente opacos
                        Screen('BlendFunction', auxWin, GL_ONE, GL_ZERO);
                        Screen('DrawTextures', auxWin, noiseTex, srcRects, dstRects, orientation(:, idx, b), [], [], [], []);
                    
            % (i) Desenha os gratings, somando ambos os sinais
                        Screen('BlendFunction', auxWin, GL_ONE, GL_ONE);
                        Screen('DrawTextures', auxWin, txP.gabor.tex, [], dstRects, orientation(:, idx, b), [], [], [], [], [], txP.gabor.props);
                    
            % (j) Desenha a abertura gaussiana
                        Screen('BlendFunction', auxWin, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                        Screen('DrawTextures', auxWin, txP.blob.tex, [], dstRects, orientation(:, idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
%                         Screen('Close', noiseTex);
            
        
            % (k) Atualiza a tela para exibir os estímulos
                        if mode == 1, Screen('DrawTexture', dpP.window, bg); lastPos = [-1 -1]; end
                        trialOnset = Screen('Flip', dpP.window, vbl + 0.5 * dpP.ifi);
            
            % vi. Registra o momento de início do trial
                        if debug == 0  && mode >= 2 && keepGoingTrials
                            Eyelink('Message',prm.msg.on.P2);
                        end
                    else 
                        trialOnset = GetSecs;
                    end
        
        % vii. O sujeito deve visitar exatamente Ns(i) estímulos
        %      antes de ocorrerem as modificações pré-sacádicas,
        %      desde que tais modificações não sejam interrompidas
        %      por movimento ocular -- nesse caso, é aumentado o
        %      valor de Ns(i)
                    counter = 0;                 % Número de estímulos visitados
                    flag = zeros(1, nStims); % Quantas vezes cada estímulo foi visitado
                    currStim = 0;
                    fixStartTime = NaN;
                    if keepGoingTrials && ~restartTrial
%                         checkFixOnset = trialOnset;
                        if debug > 0 && mode > 1
                            seenIdx = sort(randsample(nStims,modTimes(b, idx)))';
                            flag(seenIdx) = 1;
                            currIdx = randsample(seenIdx, 1);
                            KbPressWait;
                        else
                            % O loop da fase 2 só é interrompido se tiver que reiniciar o trial (por
                            % excesso de tempo) ou se o período pré-ruído rosa tiver sido completado
                            % sem movimento ocular
                            checkUpdate = false; keepP3 = true; auxModTimes = modTimes(b, idx);
                            while ~restartTrial
                                WaitSecs(0.001);
                                tNow = GetSecs;
                            
                                if checkUpdate
                                    if tNow >= preUpdateDeadline
                                        disp('De fato era último estímulo');
                                        break;
                                    end
                                end
                            
                                if auxModTimes == nStims
                                    disp('Tenho que descartar esse trial, pois excedeu o máximo de estímulos')
                                    keepP3 = false;
                                    break;
                                end
                            
                            
                                check = false;
                                if mode > 1
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
                                        Screen('Flip', dpP.window);
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
                                        if checkUpdate
                                            disp('Não era último, vamos incrementar modTimes')
                                            checkUpdate = false;
                                            auxModTimes = auxModTimes + 1;
                                        end
                                        if currStim ~= 0
                                            fixDur = tNow - fixStartTime;
                                
                                            if fixDur >= prm.minFixTime2
                                                if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.stm{1}); end
                                                % Apenas salva na fila a primeira
                                                % fixação em um estímulo
                                                if flag(currStim) == 0
                                                    counter = counter + 1;
                                                    fprintf('Terminou a visita ao %d-ésimo estímulo\n', counter)
                                                    auxFixQueue(counter) = fixDur;
                                                    P3On = P3Onset2(tkP, prm, fixDur);
                                                end
                                                seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currStim; fixDur]]; % Se quisesse registrar o comprimento de todas as fixações
                                                flag(currStim) = flag(currStim) + 1;
                                                % Registra como ruim a fixação se tiver sido muito curta
                                            else
                                                if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.stm{2}); end
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
                            
                                            %% IMPORTANTE: Se for iniciada a modTimes(b,i)-ésima fixação
                                            % diferente num estímulo, começa a contar o tempo de
                                            % atualização
                                            if flag(currStim) == 0 && counter == auxModTimes - 1
                                                checkUpdate = true;
                                                preUpdateDeadline = fixStartTime + P3On;
                                                fprintf('Visita ao (potencialmente) último %d-ésimo estímulo\n', counter+1)
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
                                    
                                    Screen('Flip', dpP.window);
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
                            end
                        end
                    end
                    % Foram mais para baixo
%                     seenIdx = find(flag ~= 0);
%                     notSeenIdx = find(flag == 0);
                    fprintf('currStim final = %d\n', currStim);
        
    
    %% 7) Início Fase 3: tela com ruído rosa
        % Apenas se a tela
        
        % (l) Desenha ruído orientado em todos os estímulos menos o 
        %     atual -- as linhas comentadas servem para mudar apenas os 
        %     não vistos. Atrasa a apresentação para o estímulo durar
        %     medFixTime segundo antes de o ruído rosa substituí-lo
                    if mode == 1
                        bg1 = Screen('OpenOffscreenWindow', dpP.window, drP.grey);
                        auxWin  = bg1;
                    else
                        auxWin = dpP.window;
                    end
                    if ~restartTrial && keepGoingTrials
                        if keepP3
                            blinkIdx = setdiff(1:(nStims), currIdx);
    
                % (i) Desenha os gratings, somando ambos os sinais
                            Screen('BlendFunction', auxWin, GL_ONE, GL_ZERO);
                %             Screen('DrawTextures', dpP.window, oriPinkTex, srcRects(:,notSeenIdx), dstRects(:,notSeenIdx), orientation(notSeenIdx, idx, b), [], [], [], []);
                            Screen('DrawTextures', auxWin, oriPinkTex, srcRects(:,blinkIdx), dstRects(:,blinkIdx), orientation(blinkIdx, idx, b), [], [], [], []);
    
                            if ~isempty(currIdx)
                                Screen('DrawTextures', auxWin, noiseTex, srcRects(:, currIdx), dstRects(:, currIdx), orientation(currIdx, idx, b), [], [], [], []);
                                Screen('BlendFunction', auxWin, GL_ONE, GL_ONE);
                                Screen('DrawTexture', auxWin, txP.gabor.tex, [], dstRects(:, currIdx), orientation(currIdx, idx, b), [], [], [], [], [], txP.gabor.props);
                            end
                        
                % (j) Desenha a abertura gaussiana
                            Screen('BlendFunction', auxWin, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                %             Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects(:,notSeenIdx), orientation(notSeenIdx, idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
                            Screen('DrawTextures', auxWin, txP.blob.tex, [], dstRects, orientation(:, idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
                            Screen('Close', oriPinkTex); Screen('Close', noiseTex);
                            Screen('BlendFunction', auxWin, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                            
                            %% A tela não modificada deve ser exibida ao todo por 
                            % medFixTime - prm.pinkNoiseDur, mas devo descontar
                            % o tempo passado desde o início da fixação
    %                         P3On = (medFixTime - prm.pinkNoiseDur);#
    %                         auxT2 = (GetSecs - fixStartTime);
    %                         timeLeft = max(0, preUpdateDur - auxT2);
                            fprintf('Tempo permitido de fixação antes do rosa: %.4f\n', P3On)
    %                         fprintf('Tempo transcorrido desde início da fixação: %.4f\n', auxT2)
    %                         fprintf('Tempo restante ate exibir ruído rosa: %.4f\n', timeLeft)
                    
                % viii. Registra os tempos de início e fim da Fase 3
    %                         preUpdateDeadline = fixStartTime + P3On;
    %                         if mode == 1
    %                             lastPos = [-1 -1];
    %                             while GetSecs < preUpdateDeadline
    %                                 [x_gaze, y_gaze, ~] = GetMouse(dpP.window);
    %                                 if any([x_gaze, y_gaze] ~= lastPos)
    %                                     Screen('DrawTexture', dpP.window, bg);
    %                                     Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
    %                                     Screen('Flip', dpP.window);
    %                                     lastPos = [x_gaze, y_gaze];
    %                                 end
    %                                 WaitSecs(0.001);
    %                             end
    % 
    %                             updateStimOnset = GetSecs;
    %                         else
    %                             updateStimOnset = Screen('Flip', dpP.window, preUpdateDeadline);
    %                         end
                            if mode == 1
                                Screen('DrawTexture', dpP.window, bg);
                                Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
                            end
                            
                            updateStimOnset = Screen('Flip', dpP.window);
                            if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.on.P3); end
    
    
                            if mode == 1
                                Screen('Close', bg);
                                clear auxWin bg;
                                % img = Screen('GetImage', dpP.window); bg = Screen('MakeTexture', dpP.window, img);
                            end
                            
                            % A fase 3 é encerrada se o estímulo fica tempo
                            % demais na tela ou quando o olho sai do último
                            % estímulo
                            fprintf('currStim último fixado (agora em P3) = %d\n', currStim);
                            while true 
                                tNow = GetSecs;
                                if tNow - updateStimOnset > prm.pinkNoiseDur
                                    P3Dur = tNow - updateStimOnset;
                                    fixDur = tNow - fixStartTime;
                                    fprintf('Fim ruído rosa por duração: %.4f\n', P3Dur)
                                    if debug == 0 && mode >= 2,  Eyelink('Message',prm.msg.off.stm{3}); end
                                    break;
                                end
                                check = false;
                                if mode > 1
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
                                        Screen('DrawTexture', dpP.window, bg1);
                                        Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
                                        Screen('Flip', dpP.window);
                                        lastPos = [x_gaze, y_gaze];
                                    end
                                end
                                if check
                                    isCurrStim = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, idx, b)) <= minFixDist1;
                                    if isCurrStim(currStim) == 0
                                        P3Dur = tNow - updateStimOnset;
                                        fixDur = tNow - fixStartTime;
                                        fprintf('Fim ruído rosa por dispersão: %.4f\n', P3Dur)
                                        if debug == 0 && mode >= 2,  Eyelink('Message',prm.msg.off.stm{1}); end
                                        break;
                                    end
                                end
                                WaitSecs(.0005);
                            end
                            % Como saio do loop da fase 2 assim que inicio a
                            % última fixação, tenho que adicionar neste momento
                            % a fixação iniciada lá
                            seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currStim; fixDur]];
                            flag(currStim) = flag(currStim) + 1;
    
                            Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                            updateStimOffset = Screen('Flip', dpP.window);
                            if debug == 0 && mode >= 2,  Eyelink('Message',prm.msg.off.P3); end
    
                            P3Dur = updateStimOffset - updateStimOnset;
                            fprintf('Tempo total de ruído rosa: %.4f\n', P3Dur)
                        end

                        seenIdx = find(flag ~= 0);
                        notSeenIdx = find(flag == 0);
                        
%                         updateStimOffset = Screen('Flip', dpP.window, updateStimOnset + prm.pinkNoiseDur);
            
            %% ix. Verifica se há alguma fixação de duração mínima em estímulo 
            %     numa janela pós-modificação
                        if keepP3
                            fixOnset = updateStimOffset; currIdx = [];
    
                            if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.on.PM); end
                
                            if debug == 0 || mode == 1
                                maxDurReached = true;
                                tNow = GetSecs;
                                while tNow - updateStimOffset < prm.postModDur
                                    check = false;
                                    if mode > 1
                                        if tNow - updateStimOffset > prm.blobPMDur
                                            Screen('Flip', dpP.window);
                                        end
                                            
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
                                                if fixDur >= prm.minFixTime3
                                                    seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currIdx; fixDur]];
                                                    disp(['Trial ' num2str(idx) ': Visitou o alvo ' num2str(currIdx) ' pós-modificação']);
                                                    if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.stm{1}); end
                                                    maxDurReached = false;
                                                    break
                                                end
                                                % Se a fixação não foi longa o suficiente, pode ser
                                                % que estava apenas passando pelo estímulo para 
                                                % fovear outro, então recomeça a busca
                                                currIdx = [];
                                                if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.stm{2}); end
                                            end
                                        end
    
    
                                        if mode == 1
                                            [x_gaze, y_gaze, ~] = GetMouse(dpP.window);
                                            if any([x_gaze, y_gaze] ~= lastPos)
                                                if tNow - updateStimOffset <= prm.blobPMDur
                                                    Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation(:, idx, b), [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                                                end
                                                Screen('FillOval', dpP.window, drP.white, [x_gaze-prm.cursorRadius_px y_gaze-prm.cursorRadius_px x_gaze+prm.cursorRadius_px y_gaze+prm.cursorRadius_px]);
                                                Screen('Flip', dpP.window);
                                                lastPos = [x_gaze, y_gaze];
                                            end
                                        end
                                    end
                                    WaitSecs(.0005);
                                    tNow = GetSecs;
                                end
                                if maxDurReached && ~isempty(currIdx)
                                   seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currIdx; prm.postModDur]];
                                   disp(['Trial ' num2str(idx) ': Visitou o alvo ' num2str(currIdx) ' pós-modificação (fim forçado)']);
                                   if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.stm{4}); end
                                end
                            end
                            fprintf('seenStimsQueue final: '); disp(seenStimsQueue{b,i}(1,:)); fprintf('\n')
    
                            if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.PM); end
                            Screen('Flip', dpP.window);
                        end
                
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
                            Screen('Close', bg1); clear auxWin bg1;
                        end

                        % Se houver fixado em algum estímulo em PM...
                        if ~isempty(currIdx)
                            nPre = nStimsToReport(1, idx, b); nPost = nStimsToReport(3, idx, b);
                        
                        % ... garante que serão perguntadas as orientações
                        % de 2 ou 3 estímulos sempre (nunca só 1)
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
                            nPre = nStimsToReport(1, idx, b); nPost = nStimsToReport(3, idx, b);


                        % ... e substitui um dos a ser perguntado
                        % (preferencialmente o já visto)
                            fprintf('Há currIdx... ')
                            nStimsToReport(2, idx, b) = 1;
                            if nPre+nPost == 2
                                nStimsToReport(1, idx, b) = max(0, nStimsToReport(1, idx, b) - 1);
                            elseif nPre+nPost >=3
                                if nPre > nPost, nStimsToReport(1, idx, b) = nPre-1;
                                else, nStimsToReport(3, idx, b) = nPost-1;
                                end
                            end
                            fprintf('então será perguntado sobre\n');
                            disp(nStimsToReport(:, idx, b))
                        end

                        allTargets = nan(1,nStims); allColors2 = drP.allColors;
                        Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                        % Divide em vizinhanças os demais pontos que não o
                        % fixado antes da atual
                        [nbhd1, nbhd2, nbhdElse] = getHexNeighborhoods2(stimCenters(:, :, idx, b), currStim, rMax);

                        % O do passado não importa de onde eu pergunto
                        auxIdx = setdiff(seenIdx, currIdx);
                        seenAux = datasample(auxIdx, min(length(auxIdx), nStimsToReport(1, idx, b)), 'Replace', false);
                        currAux = []; 
                        
                        if nStimsToReport(2, idx, b) == 1, currAux = currIdx; end

                        % O não visto eu pergunto preferencialmente da
                        % primeira vizinhança
                        chooseNbhd1 = prm.nbhd1Perc >= rand;
                        if chooseNbhd1
                            auxNotSeenIdx = intersect(notSeenIdx, nbhd1); 
                            nSnbhd(b, i) = 1;
                            disp('Vai perguntar da vizinhança próxima')
                        else
                            auxNotSeenIdx = intersect(notSeenIdx, nbhd2);
                            nSnbhd(b, i) = 2;
                            disp('Vai perguntar da vizinhança segunda')
                        end

                        if isempty(auxNotSeenIdx)
                            auxNotSeenIdx = intersect(notSeenIdx, nbhdElse);
                            nSnbhd(b, i) = 3;
                            disp('Vai perguntar da vizinhança distante')
                        end
                        
                        if currStim ~= 0
                            fprintf('Distâncias (rMax = %.4f): ', rMax);
                            disp(vecnorm(stimCenters(:, currStim, idx, b) - stimCenters(:, auxNotSeenIdx, idx, b)));
                        end

                        if isempty(auxNotSeenIdx)
                            nSnbhd(b, i) = 0;
                        end

                        notSeenAux = datasample(auxNotSeenIdx, min(length(auxNotSeenIdx), nStimsToReport(3, idx, b)), 'Replace', false);
                        % Tanto seenAux como notSeenAux devem ser linhas
                        % para concatenar 
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
                        
                                allTargets(orderToReportStims(j)) = currTarget;
                        
                                foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, rectColors, allTargets, targetOri(b), rectPW);
                        
                            end
                        end
                        
                        if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.P4); end

                        feedback = (orientation(:,idx,b) == targetOri(b))' + allTargets;
                        feedback(rem(feedback,2) == 0) = 2; feedback(rem(feedback,2) == 1) = 0;
                        feedback = feedback/2;

                        trialFeedback{b, i} = [orderToReportStims; orderRemapped; feedback(orderToReportStims)];
                        
                        if mode < 3
                            foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW, feedback, drP.red, drP.green);
                        else
                            foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW);
                         end
                        WaitSecs(.5);
                        
                        tkP.fixQueue = [tkP.fixQueue(auxModTimes:end), auxFixQueue(1:(auxModTimes-1))];
                        if keepP3
                            tkP.fixProps.preP3 = [tkP.fixProps.preP3 P3On];
                            tkP.fixProps.med   = [tkP.fixProps.med    med];
                            trialOrder(2, i, b) = 1;
    
                            i = i + 1;
                            trialIdxUp = true;
                            modTimes(b, idx) = auxModTimes;
                        else
                            if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.err.trl{3}); end
                        end

                    % A condição é keepGoingTrials = false ou restartTrials = true
                    end

                    if ~(~restartTrial && keepGoingTrials) || ~keepP3
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
                        trialOffset = GetSecs; %#ok<NASGU>
                        Eyelink('Message', sprintf(prm.msg.off.trl{1}, i - trialIdxUp, tkP.nTrials));
                        Screen('Flip', dpP.window);
                        WaitSecs(0.1);
                        Eyelink('SetOfflineMode');
                        Eyelink('StopRecording');
                    else 
                        Screen('Flip', dpP.window);
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

            if mode > 1
                tkS(mode - 1, 2) = blocksCompleted;
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
            results.trialFeedback = trialFeedback;
            results.seenStimsQueue = seenStimsQueue;
            
            if debug ~= 0, tkS(:) = 0; end
            foragingSave(tkS, 2, prm, dpP, drP, tkP, txP, results);
            cleanup(dpP.window);
            if debug == 0, diary off; end
            if sum(tkS(:)) == 0 && exist(tkP.winOutFile, 'file')
                disp('Deletando registros em texto...')
                delete(tkP.winOutFile);
            end
            psychrethrow(psychlasterror);
        end

    if mode < 3
        tkP.nBlocks = nBlocks;
        tkP.nTrials = nTrials;
        if debug == 0, HideCursor(dpP.window); end

        % Se estiver no modo treino, ajusta a duração do ruído rosa
        % conforme as fixações salvas na fila
        if mode == 2
            disp('Ajuste de pinkNoiseDur')
            shortFixLimit = prctile(tkP.fixQueue, prm.shortFixPerc);
            tkP.pinkNoiseDur = min(prm.pinkNoiseDur, shortFixLimit);
        end
    end
end



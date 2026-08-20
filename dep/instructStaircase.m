function instructStaircase(tkP, dpP, drP, txP, prm)
% Instrução mostra telas principais e algumas adicionais:
% 1. Mostrando que a qualidade dos estímulos pode mudar entre trials
% 2. Mensagens de erro

    %% Inicializa variáveis úteis
    Screen('Flip', dpP.window);
    nBlocks = 1; nTrials = 1;
    
    leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4}; sKey = tkP.keys{6};
    if strcmp(tkP.targetKey, 'right'), targetKey = rightKey; % nonTargetKey = leftKey;
    else,                              targetKey = leftKey;  % nonTargetKey = rightKey;
    end
    
    nStims = prm.nStimsStair;
    [nTs,  ~, targetOri, ~, ~, ~] = getForagingDistributions1(nStims, tkP.nMinFix, tkP.nMaxFix, nTrials, nBlocks, prm);
    
    drP.allColors = drP.white*ones(3, nStims); drP.allPW     = prm.pW1*ones(1,nStims);
    
    lowerBound = targetOri + prm.nbhdRadius;
    upperBound = 180 + (targetOri - prm.nbhdRadius);
    orientation = mod(rand(nStims, nTrials)*(upperBound-lowerBound)+lowerBound, 180);

    crossSize_px = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.crossSize_dva);

    noiseMatrix = zeros(txP.gabor.size_px, (nStims)*txP.gabor.size_px);
    oriPinkMatrix = zeros(txP.gabor.size_px, (nStims)*txP.gabor.size_px);

    noiseCenters = [0:txP.gabor.size_px:(nStims-1)*txP.gabor.size_px; zeros(1, nStims)] + txP.gabor.size_px/2;
    baseRect = [0 0 txP.gabor.size_px txP.gabor.size_px];

    srcRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', noiseCenters(1,:), noiseCenters(2,:));

    minDist_px   = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.minDist_dva+prm.gaborSize_dva);

    [currFixCenter, currStimCenter, ~, fI] = getStimLocations2_1(dpP.winRect(3:4), [dpP.winCenter 1], tkP.nStims, minDist_px, txP.gabor.size_px);
    fixCenters     = currFixCenter;
    aux = NeighborsOrder(currStimCenter, fI, []);
    stimCenters = currStimCenter(:, [aux(1:nStims-1), fI]);
    fixIdx = nStims;
    orientation(randperm(nStims, nTs)) = targetOri;

    stimRange = (-prm.sigmaMax:.5:-prm.sigmaMin);
    [~, auxidx] = min(abs(stimRange - (-prm.aSigma(1))));
    prm.aSigma(1) = stimRange(auxidx);

    [oriFilter, OFsize] = MakeOriFilter1(txP.gabor.size_px, prm.aSigma(1), prm.rSigma2);

    b = 1;
    dstRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', stimCenters(1,:), stimCenters(2,:));
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


%     auxNoiseMatrix = butterFilter(pinkNoise(txP.gabor.size_px, (nStims)*txP.gabor.size_px), txP.noiseLoCutFreq, txP.noiseHiCutFreq);
    auxNoiseMatrix = pinkNoise(txP.gabor.size_px, (nStims)*txP.gabor.size_px);
    for j=1:(nStims)
        colRange = ((j-1)*txP.gabor.size_px+1):(j*txP.gabor.size_px);
        aux = auxNoiseMatrix(:,colRange);
        aux1 = butterFilter(aux, txP.noiseLoCutFreq, txP.noiseHiCutFreq);
        noiseMatrix(:, colRange) = (aux1 - mean(aux1(:)))/std(aux1(:));
        oriPinkMatrix(:,colRange) = ApplyOriFilter1(oriFilter', OFsize, aux);
    end
    noiseTex   = Screen('MakeTexture', dpP.window,  prm.noiseSTDmult*noiseMatrix,      [], [], 1);
    gaborTex   = Screen('MakeTexture', dpP.window,  prm.gaborSTDmult*txP.gabor.matrix, [], [], 1);
    oriPinkTex = Screen('MakeTexture', dpP.window,  prm.stimSTDmult *oriPinkMatrix,    [], [], 1);

    xFix = [-crossSize_px/2 crossSize_px/2 0 0];
    yFix = [0 0 -crossSize_px/2 crossSize_px/2];
    fixCoords = [xFix; yFix];

    blinkIdx = setdiff(1:nStims, fixIdx);
    alphas = ones(nStims, 1);

    nbhd = NeighborsOrder(stimCenters, fixIdx, []);
    stimsToReport = nbhd(1:nStims-1);
    orderToReportStims = stimsToReport(randperm(nStims-1));

    allTargets = nan(1,nStims);

    overlayRect = CenterRectOnPointd([dpP.winRect(1:3) dpP.winRect(4)/2.5], dpP.winRect(3)/2, dpP.winRect(4)/2);
    overlayColor = repmat(drP.whiteGrey, 1, 3);
    textColor1   = repmat(drP.black, 1, 3);
    textColor2   = repmat(.5*(drP.grey+drP.black), 1, 3);

    eyeCols = [drP.blue; drP.darkBlue; drP.darkGreen; drP.brown; drP.darkBrown; drP.greyBrown; drP.paleBrown]; aux = randsample(1:size(eyeCols, 1), 1);
    eyeColor = eyeCols(aux, :);

    interestOri  = 90:45:225;
    interestFilt = 5 * 2.^(0:4);
    allPinkTex = createNewPinkTex(dpP, txP, prm, interestOri, interestFilt);

    HideCursor;
    selected = 1; 
%     fromAtypical = false;
%     prevState = 0;
%     atypicalIdx = [11; 12];
%     gotoAtypical = [3 5];
    %% Tarefa
    
    screens = { 
        @(x)drawBlocks(tkP, dpP, drP, txP, prm, targetOri, b);                          % Tela 1: blocos 
        @(x)DrawCalibrationTarget(dpP.window, dpP.winCenter(1), dpP.winCenter(2));      % Tela 2: calibração
        @(x)startFix(tkP, dpP, drP, txP, prm, fixCoords, fixCenters);                   % Tela 3: fixação inicial
        @(x)foragingDrawMain(dpP.window, gaborTex, noiseTex, srcRects, ...
        dstRects, orientation, txP, [repmat(alphas', [3,1]); ones(1, nStims)]);         % Tela 4: estímulos
        @(x)drawPink(tkP, dpP, drP, txP, prm, oriPinkTex, gaborTex, ...
        noiseTex, srcRects, dstRects, blinkIdx, fixIdx, orientation);                   % Tela 5: ruído rosa
        @(x)drawPM(tkP, dpP, drP, txP, prm, dstRects, orientation, textColor2);         % Tela 6: pedestais
        @(x)drawInteractive(tkP, dpP, drP, txP, prm, nStims, orderToReportStims, ...    % Tela 7: interativa de resposta
        leftKey, rightKey, spaceKey, dstRects, orientation, textColor2, stimCenters, targetOri);
        @(x)drawPostAns(tkP, dpP, drP, txP, prm, dstRects, orientation, ...
        textColor2, stimCenters, orderToReportStims, x, targetOri);                     % Tela 8: pós-resposta
        @(x)obsPink(tkP, dpP, drP, txP, prm, x, gaborTex, ...                           % Tela 9: diferenças de intensidade
        noiseTex, interestOri)
        @(x)obsArrows(tkP, dpP, drP, txP, prm, leftKey, rightKey, targetKey, ...        % Tela 10: correspondência entre 
        targetOri)                                                                      %          setas e alvo ou distrator
        @(x)startFix1(tkP, dpP, drP, txP, prm, fixCoords, fixCenters);                  % Tela 11
        @(x)warnTimeOut(tkP, dpP, drP, txP, prm, fixCoords, fixCenters, ...             % Tela 12: aviso de tempo  (antecedida por tela 3/11)
        overlayColor, overlayRect, textColor1)  ;
        @(x)drawPink(tkP, dpP, drP, txP, prm, oriPinkTex, gaborTex, ...
        noiseTex, srcRects, dstRects, blinkIdx, fixIdx, orientation);                   % Tela 13: aviso de que mexeu o olho (antecedida por tela 5/13) 
        @(x)warnMov(tkP, dpP, drP, txP, prm, stimCenters, dstRects, fixIdx);
        };

    nScreens = numel(screens);

    % Posiciona o cursor no centro mas ligeiramente deslocado
    r = 5*prm.eyeSize;
    theta = 2*pi*rand; dx = r*cos(theta); dy = r*sin(theta);
    SetMouse(dpP.winCenter(1) + dx, dpP.winCenter(2) + dy, dpP.window);

    oldState = false;
    oldSState = false;
    quit = false;

    [~,~, buttons] = GetMouse(dpP.window);

    while ~quit
        Screen('FillRect', dpP.window, drP.grey);

        % Desenha a tela selecionada

        if selected == 7
            [outKey, allTargets] = screens{selected}(dpP.window);
            pressed = true; 
            kc = zeros(1, numel(kc));kc(outKey) = 1;
            oldState = 0;
        else
            if selected == 9
                screens{selected}(allPinkTex);
            else
                screens{selected}(allTargets);
            end

            [mx,my, buttons] = GetMouse(dpP.window);
            DrawEye(dpP.window, [mx my], prm.eyeSize, eyeColor);
    
            Screen('Flip', dpP.window);

            [pressed,~,kc] = KbCheck;

            if selected == 9 && kc(sKey) && ~oldSState
                allPinkTex = createNewPinkTex(dpP, txP, prm, interestOri, interestFilt, allPinkTex);
            end
            oldSState = kc(sKey);
        end

        if (pressed || any(buttons)) && ~oldState
            if kc(escapeKey)
                quit = true;

            elseif kc(rightKey) || (length(buttons) >= 3 && buttons(3))
                if selected < nScreens
                    selected = selected + 1;
                else
                    % Encerra as instruções se for a última tela
                    quit = true;
                end
            elseif kc(leftKey) || buttons(1)
                if selected == 12 || selected == 14, selected = selected - 2;
                else,                                selected = max(1, selected-1);
                end
            end
        end
        oldState = (pressed || any(buttons));

        WaitSecs('YieldSecs', 0.005);
    end
    KbWait([], 1);
    Screen('Close', noiseTex);    Screen('Close', gaborTex); Screen('Close', oriPinkTex)
    for i = 1:numel(interestOri)
        for j = 1:numel(interestFilt)
            Screen('Close', allPinkTex(i,j));
        end
    end
end

function drawBlocks(tkP, dpP, drP, txP, prm, tgtOri, b)
    Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
    Screen('DrawTexture', dpP.window, txP.exampleGabor.tex, [], [], tgtOri, [], [], [], [], []);
    Screen('DrawTexture', dpP.window, txP.exampleNoise.tex, [], [], [], [], [], [], [], []);
    Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
    Screen('DrawTextures', dpP.window, txP.exampleBlob.tex, [], [], tgtOri, [], [], [0 0 0 1]', [], [], txP.exampleBlob.props); 

    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('TextSize', dpP.window, prm.textSizeBigger);
    countBlockText  = sprintf('Bloco %d de %d', b, tkP.nBlocks);
    DrawFormattedText(dpP.window, countBlockText, 'center', dpP.screenRes.height*.1, drP.black);

    Screen('TextSize', dpP.window, prm.textSizeHuge);
    targetText_1    = 'Orientação dos alvos:';
    DrawFormattedText(dpP.window, targetText_1, 'center', dpP.screenRes.height*.2, drP.black);

    targetText_2    = prm.allOriName{prm.allOriMap(tgtOri)};
    Screen('TextStyle', dpP.window, 1); Screen('TextSize', dpP.window, prm.textSizeHuger);
    DrawFormattedText(dpP.window, targetText_2, 'center', dpP.screenRes.height*.3, drP.black);
                                                    
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('TextStyle', dpP.window, 0);
    Screen('TextSize', dpP.window, prm.textSizeHuge);
    proceedText = 'Pressione ESPAÇO para prosseguir';
    DrawFormattedText(dpP.window, proceedText, 'center', dpP.screenRes.height*.75, drP.black);

    Screen('TextSize', dpP.window, prm.textSizeNormal);

    Screen('TextFont', dpP.window, prm.textFont);
end

function startFix(~, dpP, drP, ~, prm, fixCoords, fixCenters)
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('DrawLines', dpP.window, fixCoords, prm.lineWidth_px, drP.white, fixCenters', 2);
end

function startFix1(~, dpP, drP, ~, prm, fixCoords, fixCenters)
    Screen('TextSize', dpP.window, prm.textSizeEnormous);
    titleText = 'Advertências'; titleMargin = dpP.winCenter(1)*prm.titleMarginFactor;
    Screen('DrawText', dpP.window, titleText, titleMargin, titleMargin,  drP.blackGrey);

    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('DrawLines', dpP.window, fixCoords, prm.lineWidth_px, drP.white, fixCenters', 2);

    Screen('TextSize', dpP.window, prm.textSizeNormal);
    Screen('TextFont', dpP.window, prm.textFont);
end

function drawPink(~, dpP, ~, txP, ~, oriPinkTex, gaborTex, noiseTex, srcRects, dstRects, blinkIdx, fixIdx, orientation)
    Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
    if isempty(blinkIdx)
        Screen('DrawTextures', dpP.window, oriPinkTex, srcRects, dstRects, orientation, [], [], [], []);
    else
        Screen('DrawTextures', dpP.window, oriPinkTex, srcRects(:,blinkIdx), dstRects(:,blinkIdx), orientation(blinkIdx), [], [], [], []);
    end
    
    if ~isempty(fixIdx)
        Screen('DrawTextures', dpP.window, gaborTex, [], dstRects(:,fixIdx), orientation(fixIdx));
        Screen('DrawTextures', dpP.window, noiseTex, srcRects(:,fixIdx), dstRects(:,fixIdx), orientation(fixIdx));
        
        Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
        Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects(:, fixIdx), orientation(fixIdx), [], [], [0 0 0 1]', [], [], txP.blob.props);
    end
    Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
    if isempty(blinkIdx)
        Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects, orientation, [], [], [0 0 0 1]', [], [], txP.blob.props);
    else
        Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects(:, blinkIdx), orientation(blinkIdx), [], [], [0 0 0 1]', [], [], txP.blob.props);
    end
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
end

function drawPM(~, dpP, ~, txP, ~, dstRects, orientation, textColor2)
    Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation, [], [], [textColor2 1]', [], [], txP.PMBlob.props);
end

function [outKey, allTargets] = drawInteractive(tkP, dpP, drP, txP, prm, nStims, orderToReportStims, leftKey, rightKey, spaceKey, dstRects, orientation, textColor2, stimCenters, targetOri)
    if strcmp(tkP.targetKey, 'right'), targetKey = rightKey; nonTargetKey = leftKey;
    else,                              targetKey = leftKey; nonTargetKey = rightKey;
    end
    allTargets = nan(1,nStims); allColors2 = drP.allColors;
    N = numel(orderToReportStims);

    j=1;
    rectColors = allColors2; rectColors(:, orderToReportStims(j)) = drP.orange;rectPW = drP.allPW; rectPW(:, orderToReportStims(j)) = prm.pW2;

    Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation, [], [], [textColor2 1]', [], [], txP.PMBlob.props);
    foragingFlip(dpP.window, stimCenters, dstRects, orderToReportStims, txP.gabor.size_px, rectColors, allTargets, targetOri, rectPW);
    KbReleaseWait;
    
    while true
        [~,~,buttons] = GetMouse(dpP.window);
        if ~any(buttons)
            break;
        end
        WaitSecs(0.001);
    end
    keyPressed = zeros(1, max([rightKey, leftKey, spaceKey]));
    [~,~, buttons] = GetMouse(dpP.window);
    buttons(:) = 0;

    while ~(any(keyPressed([spaceKey, rightKey, leftKey])) || any(buttons))
        while true
            [~,~,keyPressed] = KbCheck;
            [~, ~, buttons] = GetMouse(dpP.window);
            if any(keyPressed) || any(buttons)
                KbWait([], 1);
                break;
            end
            WaitSecs(0.001);
        end
        if ~any(buttons) && keyPressed(spaceKey)
            outKey = rightKey;
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
                    Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation, [], [], [textColor2 1]', [], [], txP.PMBlob.props);
                    foragingFlip(dpP.window, stimCenters, dstRects, orderToReportStims, txP.gabor.size_px, rectColors, allTargets, targetOri, rectPW);
                end
            end
        elseif keyPressed(rightKey) || (length(buttons) >= 3 && buttons(3))
            outKey = rightKey;
            allTargets(orderToReportStims) = randi(2, [1 N])-1;
        elseif keyPressed(leftKey) || buttons(1)
            outKey = leftKey;
            allTargets(orderToReportStims) = randi(2, [1 N])-1;
        end
        WaitSecs('YieldSecs', 0.005);
    end
    KbReleaseWait;
end

function drawPostAns(~, dpP, drP, txP, ~, dstRects, orientation, textColor2, stimCenters, orderToReportStims, allTargets, targetOri)

    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('DrawTextures', dpP.window, txP.PMBlob.tex, [], dstRects, orientation, [], [], [textColor2 1]', [], [], txP.PMBlob.props);
    foragingFlip1(dpP.window, stimCenters, dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri, drP.allPW);
end

function obsPink(tkP, dpP, drP, txP, prm, allPinkTex, gaborTex, noiseTex, interestOri)

    [w, h] = RectSize(dpP.winRect);

    % Grid dimensions
    [nRows, nCols] = size(allPinkTex);

    % Image dimensions
    imageWidth  = txP.gabor.size_px;
    imageHeight = txP.gabor.size_px;

    % Spacing between images
    xSpacing = imageWidth  * 2.5;
    ySpacing = imageHeight * 1.8;

    % Center the grid on the screen
    gridCenterX = w / 2;
    gridCenterY = 4*h / 7;

    % Total grid dimensions
    gridWidth  = (nCols - 1) * xSpacing + imageWidth;
    gridHeight = (nRows - 1) * ySpacing + imageHeight;

    % Top-left position of the grid
    gridX0 = gridCenterX - gridWidth / 2;
    gridY0 = gridCenterY - gridHeight / 2;

    % Draw grid
    for row = 1:nRows
        for col = 1:nCols

            x = gridX0 + (col - 1) * xSpacing;
            y = gridY0 + (row - 1) * ySpacing;

            rect = CenterRectOnPointd( ...
                [0 0 imageWidth imageHeight], ...
                x + imageWidth / 2, ...
                y + imageHeight / 2);

            drawPink(tkP, dpP, drP, txP, prm, ...
                allPinkTex(row, col), gaborTex, noiseTex, [], ...
                rect, [], [], interestOri(row));

        end
    end

    % Text
    Screen('TextSize', dpP.window, prm.textSizeEnormous);

    titleText = 'Observações';
    titleMargin = dpP.winCenter(1) * prm.titleMarginFactor;

    Screen('DrawText', dpP.window, titleText, ...
        titleMargin, titleMargin, drP.blackGrey);

    Screen('TextSize', dpP.window, prm.textSizeNormalish);

    leftText  = 'Muito visível';
    rightText = 'Pouco visível';

    % X coordinates of leftmost and rightmost columns
    leftColumnX  = gridX0 + imageWidth / 2;
    rightColumnX = gridX0 + (nCols - 1) * xSpacing + imageWidth / 2;

    % Text above leftmost column
    bounds = Screen('TextBounds', dpP.window, leftText);
    textWidth = RectWidth(bounds);

    DrawFormattedText(dpP.window, leftText, ...
        leftColumnX - textWidth / 2, ...
        gridY0 - 50, ...
        drP.black);

    % Text above rightmost column
    bounds = Screen('TextBounds', dpP.window, rightText);
    textWidth = RectWidth(bounds);

    DrawFormattedText(dpP.window, rightText, ...
        rightColumnX - textWidth / 2, ...
        gridY0 - 50, ...
        drP.black);

    Screen('TextSize', dpP.window, prm.textSizeNormal);
    Screen('TextFont', dpP.window, prm.textFont);

end

function obsArrows(~, dpP, drP, txP, prm, leftKey, rightKey, targetKey, targetOri)

    [w, h] = RectSize(dpP.winRect);
    leftX  = round(w*(3/8));
    rightX = round(w*(5/8));
    
    imageWidth  = txP.gabor.size_px;
    imageHeight = txP.gabor.size_px;
    
    leftRect  = CenterRectOnPointd([0 0 imageWidth imageHeight], leftX, h/2);
    rightRect = CenterRectOnPointd([0 0 imageWidth imageHeight], rightX, h/2);

    foragingFlip1(dpP.window, [leftX, h/2]', leftRect', 1, txP.gabor.size_px, drP.orange', leftKey == targetKey, targetOri, drP.allPW(1));
    foragingFlip1(dpP.window, [rightX, h/2]', rightRect', 1, txP.gabor.size_px, drP.orange', rightKey == targetKey, targetOri, drP.allPW(1));

    Screen('TextSize', dpP.window, prm.textSizeEnormous);
    titleText = 'Observações'; titleMargin = dpP.winCenter(1)*prm.titleMarginFactor;
    leftText = 'Seta esquerda <';
    rightText = '> Seta direita';
    Screen('DrawText', dpP.window, titleText, titleMargin, titleMargin,  drP.blackGrey);
    Screen('TextSize', dpP.window, prm.textSizeNormalish);

    bounds = Screen('TextBounds', dpP.window, leftText); textWidth = RectWidth(bounds);
    DrawFormattedText(dpP.window, leftText,   leftX - textWidth/2, leftRect(4)+50, drP.black);
    bounds = Screen('TextBounds', dpP.window, rightText); textWidth = RectWidth(bounds);
    DrawFormattedText(dpP.window, rightText, rightX - textWidth/2, rightRect(4)+50, drP.black);
    
    Screen('TextSize', dpP.window, prm.textSizeNormal);
    Screen('TextFont', dpP.window, prm.textFont);
end

function warnTimeOut(~, dpP, drP, ~, prm, fixCoords, fixCenters, overlayColor, overlayRect, textColor1)
    Screen('TextSize', dpP.window, prm.textSizeEnormous);
    titleText = 'Advertências'; titleMargin = dpP.winCenter(1)*prm.titleMarginFactor;
    Screen('DrawText', dpP.window, titleText, titleMargin, titleMargin,  drP.blackGrey);

    alpha = 1;
    Screen('DrawLines', dpP.window, fixCoords, prm.lineWidth_px, drP.white, fixCenters', 2);
    
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

    Screen('TextSize', dpP.window, prm.textSizeNormal);
    Screen('TextFont', dpP.window, prm.textFont);
end

function warnMov(~, dpP, drP, txP, prm, stimCenters, dstRects, fixIdx)
    Screen('TextSize', dpP.window, prm.textSizeEnormous);
    titleText = 'Advertências'; titleMargin = dpP.winCenter(1)*prm.titleMarginFactor;
    Screen('DrawText', dpP.window, titleText, titleMargin, titleMargin,  drP.blackGrey);
    warningFlip1(dpP.window, stimCenters, resizeRect(dstRects, .25), fixIdx, txP.gabor.size_px, drP.allPW, drP.darkRed);

    Screen('TextSize', dpP.window, prm.textSizeNormal);
    Screen('TextFont', dpP.window, prm.textFont);
end


function allPinkTex = createNewPinkTex(dpP, txP, prm, interestOri, interestFilt, auxPink)
    if nargin >= 6 && ~isempty(auxPink)
        Screen('Close', auxPink(auxPink > 0));
    end

    allPinkTex = zeros(numel(interestOri), numel(interestFilt));
    
    for i = 1:numel(interestOri)
        for j = 1:numel(interestFilt)
            [easyFilter, easyOFsize] = MakeOriFilter1(txP.gabor.size_px, interestFilt(j), prm.rSigma2);
    %     auxNoiseMatrix1 = butterFilter(pinkNoise(txP.gabor.size_px, txP.gabor.size_px), txP.noiseLoCutFreq, txP.noiseHiCutFreq);
            auxNoiseMatrix1 = pinkNoise(txP.gabor.size_px, txP.gabor.size_px);
            oriPinkMatrix1  = ApplyOriFilter1(easyFilter', easyOFsize, auxNoiseMatrix1);
            allPinkTex(i,j) = Screen('MakeTexture', dpP.window,  prm.stimSTDmult *oriPinkMatrix1, [], [], 1);
        end
    end
end





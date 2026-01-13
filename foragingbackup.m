function ForagingGabors(nTrials, nStims, nBlocks, nTargets, options)
% Versão (que pretende ser) final, com estímulo pré-sacádico em alguns
% trials e report das orientações
% A tarefa consiste em nTrials trials, cada um com nStims estímulos e 
% quantidade média de nTargets alvos, espalhados pseudoaleatoriamente na 
% tela. A tarefa do sujeito é encontrar, dentre gabores de alta frequência 
% com orientações variadas, alvos com uma orientação específica, apenas com
% o movimento ocular. Em determinado momento, durante essa busca ocular, os
% estímulos ainda não visitados são substituídos por um ruído rosa de mesma 
% orientação durante o período pré-sacádico. Por fim, o sujeito deve
% reportar se viu ou não alvo nas posições indicadas
%% Depende do participante: olho registrado, orientação-alvo, fixação média
%% de cada participante, algum parâmetro do ruído rosa
    
        arguments
            nTrials {mustBeNumeric}   = 20;
            nStims {mustBeNumeric}    = 8;
            nBlocks {mustBeNumeric}   = 15;
            nTargets {mustBeNumeric}  = ceil(nStims/2);
            options.mode string       = 'experiment'
        end
        
        rng('shuffle');
        cleanup

        % Verifica se o programa será executado no modo debug ou não
        debug = 0;
        if strcmp(options.mode, 'debug'), debug = 1; elseif strcmp(options.mode, 'debugTV'), debug = 2; end
        if debug > 0
            nTrials  = 2;
            nStims   = 8;
            nBlocks  = 2;
            nTargets = ceil(nStims/2);
        end

%% 1) Importa os parâmetros para o experimento
        params = foragingParams;
        
    
%% 2) Inicializa PTB e EyelinkToolBox
    % i. Caixa de diálogo 
    if debug == 0
        prompt = {'\fontsize{15} Sujeito', '\fontsize{15} Número da sessão', '\fontsize{15} Olho dominante (E ou D)'};
        dlg_title = 'Informações da sessão experimental';
        def = {'00', '01', 'R'};
        options.Resize = 'on'; options.Interpreter = 'tex';
        answer = inputdlg(prompt, dlg_title, [1 50; 1 30; 1 30], def, options);
        domEye = answer{3};
        if  isempty(answer), fprintf('Sessão cancelada pelo usuário\n'); cleanup; return; end
        edfFile = [answer{1} '_' answer{2}];

        % Interrompe se nome muito longo
        if length(edfFile) > 8, fprintf('ERRO: O nome do arquivo excede 8 caracteres!\n'); cleanup; return; end
    end
    
    % (a) Configurações antes de abrir a tela
        if debug == 0
            FlushEvents;
            PsychDefaultSetup(2);
            Screen('Preference', 'SyncTestSettings', 0.01, 50, 0.25);
            Screen('Preference', 'SuppressAllWarnings', 1);
            Screen('Preference', 'Verbosity', 0);
            PsychImaging('PrepareConfiguration');
            PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        else
            PsychDefaultSetup(2);
            Screen('Preference', 'SkipSyncTests', 1);
        end
        KbName('UnifyKeyNames');
        leftKey   = KbName('LeftArrow'); rightKey  = KbName('RightArrow');
        spaceKey  = KbName('space');     escapeKey = KbName('ESCAPE');
        keys = {leftKey, rightKey, spaceKey, escapeKey};


    % (b) Escolhe a tela em que haverá o desenho e define algumas cores
        screenNumber = max(Screen('Screens'));
        white = WhiteIndex(screenNumber);
        grey = white / 2;
        [window, winRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2, [], [],  kPsychNeed32BPCFloat);

    % (c) Obtém propriedades da tela
        [monitorW_mm, ~] = Screen('DisplaySize', screenNumber); % Tamanho   da tela em mm
        screenRes = Screen('Resolution', screenNumber);         % Resolução da tela em px
        ifi = Screen('GetFlipInterval', window);

    % (d) Obtém o centro e os semieixos da elipse
        ellipseProps = [screenRes.width/2, screenRes.height/2];
        ellipseProps = [ellipseProps ellipseProps(1)*params.ellipseToScreenRatio(1) ellipseProps(2)*params.ellipseToScreenRatio(2)];
        Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % (e) Define outras cores de interesse e espessura das canetas para
    %     desenhar algumas das figuras
        midGrey = (white+grey)/2;
        black = BlackIndex(screenNumber);
        orange = [250, 140, 15]'/255;
        red    = [white 0 0];
        green  = [0 white 0];
        
        allColors = white*ones(3, nStims);
        allPW     = params.pW1*ones(1,nStims);

    % i. Inicializa Eyelink
        dummymode = 0; 
        if debug > 0, dummymode = 1; disp('Debug: O EyeLink não será inicializado'); end

        EyelinkInit(dummymode);

    % ii. Abre o arquivo .edf
        if debug == 0
            failOpen = Eyelink('OpenFile', edfFile);
            % Interrompe se arquivo não abrir
            if failOpen ~= 0, fprintf('Cannot create EDF file %s', edfFile); cleanup; return; end
    
    % iv. Ajusta o que é comunicado entre os PCs
            Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
            Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON,FIXUPDATE,INPUT');
    
    % v. Ajusta o que é salvo no arquivo .edf
            Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT');
            Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
    
    % vi. Identifica o olho dominante
            if ismember(domEye, {'L','E'}), Eye=1; Eyelink('Command', 'active_eye = LEFT'); elseif ismember(domEye, {'R', 'D'}), Eye=2; Eyelink('Command', 'active_eye = RIGHT'); end
    
    % vii. Adiciona um texto inicial ao arquivo .edf
            preambleText = sprintf('RECORDED BY %s session name: %s', mfilename, edfFile);
            Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);
    
    % viii. Ajusta as configurações-padrão do Eyelink de calibração
            el = EyelinkInitDefaults(window);
            el.calibrationtargetsize = 2;
            el.calibrationtargetwidth = 0.3;
            el.backgroundcolour = grey;
            el.calibrationtargetcolour = [0 0 0];
            el.msgfontcolour = [0 0 0];
            el.targetbeep = 0; el.feedbackbeep = 0;
            EyelinkUpdateDefaults(el);
            Eyelink('Command', 'screen_pixel_coords = %ld %ld %ld %ld', 0, 0, screenRes.width-1, screenRes.height-1);
            Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, screenRes.width-1, screenRes.height-1);
            Eyelink('Command', 'calibration_type = HV9'); Eyelink('command', 'generate_default_targets = NO');
            % Podemos definir o espaçamento dos pontos de calibração
            Eyelink('command', 'calibration_samples = 10'); Eyelink('command', 'calibration_sequence = 0,1,2,3,4,5,6,7,8,9');
            Eyelink('command', 'calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
                960,540, 960,205, 960,875, 442,540, 1478,540, 442,205, 1478,205, 442,875, 1478,875);
            Eyelink('command', 'validation_samples = 10'); Eyelink('command', 'validation_sequence = 0,1,2,3,4,5,6,7,8,9');
            Eyelink('command', 'validation_targets = %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
                960,540, 960,205, 960,875, 442,540, 1478,540, 442,205, 1478,205, 442,875, 1478,875);
            Eyelink('Command', 'button_function 5 "accept_target_fixation"');
    
            Eyelink('Command', 'clear_screen 0'); % Limpa a tela do Host PC
            topPriorityLevel = MaxPriority(window);
            Priority(topPriorityLevel);
            HideCursor(window);
            ListenChar(-1); % Para o teclado não escrever nada na command window
        else
            Eye = [];
        end

%% 3) Obtém as texturas de Gabor, abertura gaussiana e do ruído
    if debug == 2, params.screenDist = 200; end
    % (a) Estímulos para a tarefa
        gabor = foragingGabor(window, screenRes, nStims, params, monitorW_mm);
        blob = foragingBlob(window, gabor.size_px, grey, params);
        % Parâmetros do filtro para ruído
        noiseHiCutFreq = params.noiseHiCutFreq_cpd/gabor.ppd;
        noiseLoCutFreq = params.noiseLoCutFreq_cpd/gabor.ppd;

        % Filtro de orientação
        [oriFilter, OFsize] = MakeOriFilter(gabor.size_px, params.aSigma, params.rSigma2);

    % (b) Estímulos para a tela entre blocos
        exampleParams = params; exampleGaborFactor = 3;
        exampleParams.gaborFreq_cpd = params.gaborFreq_cpd/exampleGaborFactor;
        exampleParams.gaborSize_dva = params.gaborSize_dva*exampleGaborFactor;
        exampleGabor = foragingGabor(window, screenRes, 1, exampleParams, monitorW_mm);
        exampleBlob = foragingBlob(window, exampleGabor.size_px, grey, params);
        exampleNoise = foragingNoise(window, exampleGabor.size_px, 1, grey, params);

%% 4) Reúne todas as variáveis definidas até aqui para passá-las às demais funções

        fixQueue = nan(1, params.fixTimeQueueSize);

        displayProps.window       = window;
        displayProps.winRect      = winRect;
        displayProps.monitorW_mm  = monitorW_mm;
        displayProps.screenRes    = screenRes;
        displayProps.ifi          = ifi;
        displayProps.ellipseProps = ellipseProps;

        drawProps.white     = white;
        drawProps.grey      = grey;
        drawProps.midGrey   = midGrey;
        drawProps.black     = black;
        drawProps.orange    = orange;
        drawProps.red       = red;
        drawProps.green     = green;
        drawProps.allColors = allColors;
        drawProps.allPW     = allPW;

        texProps.gabor          = gabor;
        texProps.blob           = blob;
        texProps.exampleGabor   = exampleGabor;
        texProps.exampleBlob    = exampleBlob;
        texProps.exampleNoise   = exampleNoise;
        texProps.noiseHiCutFreq = noiseHiCutFreq;
        texProps.noiseLoCutFreq = noiseLoCutFreq;
        texProps.oriFilter      = oriFilter;
        texProps.OFsize         = OFsize;

        taskProps.Eye      = Eye;
        taskProps.fixQueue = fixQueue;
        taskProps.nTargets = nTargets;
        taskProps.nStims   = nStims;
        taskProps.nTrials  = nTrials;
        taskProps.nBlocks  = nBlocks;
        taskProps.keys     = keys;

        menuSelection(taskProps, displayProps, drawProps, texProps, debug, params);



%% 9) Desenha os Gabores em cada trial
    runForaging
    
    Eyelink('CloseFile');
    cleanup;
end

function menuSelection(tkP, dpP, drP, txP, debug, prm)


Screen('TextFont', dpP.window, 'Times');
Screen('TextSize', dpP.window, 24);

[xc, yc] = RectCenter(dpP.winRect);
SetMouse(xc, yc, dpP.window);
HideCursor;

leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4};

%% Button rectangles
btnW = 300; btnH = 100;
gap  = 80;

leftRect  = CenterRectOnPoint([0 0 btnW btnH], xc - btnW/2 - gap/2, yc);
rightRect = CenterRectOnPoint([0 0 btnW btnH], xc + btnW/2 + gap/2, yc);

% Skip arrow
skipRect = [dpP.winRect(3)-150, dpP.winRect(4)-60, dpP.winRect(3)-40, 40];

%% State
selected  = 0;
[prevMx, prevMy, prevBut] = GetMouse(dpP.window);
prevMouse = [prevMx, prevMy, prevBut];

arrow = [ ...
    -30  -8;   % tail bottom
     10  -8;   % shaft bottom
     10 -18;   % head bottom
     30   0;   % tip
     10  18;   % head top
     10   8;   % shaft top
    -30   8 ]; % tail top
scale = 4;
arrow = arrow * scale;
% Screen('FramePoly', win, color, arrow, 2);

isMouseMostRecent = false;
while true

    % --- Mouse handling ---
     [mx, my, buttons] = GetMouse(dpP.window);
    mouseMoved = norm([mx my] - prevMouse(1:2)) > 1;
    prevMouse = [mx my buttons];
    cursorVisible = false;

    if mouseMoved || isMouseMostRecent
        cursorVisible = true;
        isMouseMostRecent = true;
        if IsInRect(mx, my, leftRect)
            selected = 1;
        elseif IsInRect(mx, my, rightRect)
            selected = 2;
        else
            selected = 0;
        end
    end

    % --- Keyboard handling ---
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        if keyCode(leftKey)
            selected = 1;
            cursorVisible = false;
            isMouseMostRecent = false;
        elseif keyCode(rightKey)
            selected = 2;
            cursorVisible = false;
            isMouseMostRecent = false;
        elseif keyCode(spaceKey)
            if selected == 1
                mode = 'cursor';
            elseif selected == 2
                mode = 'train';
            end
            runForaging(tkP, dpP, drP, txP, prm, debug, mode);
            selected = 0;
        elseif keyCode(escapeKey)
            break
        end
    end

    % --- Mouse click ---
    if any(buttons)
        if IsInRect(mx, my, skipRect)
            break
        elseif selected == 1
            mode = 'cursor';
            runForaging(tkP, dpP, drP, txP, prm, debug, mode);
            selected = 0;
        elseif selected == 2
            mode = 'train';
            runForaging(tkP, dpP, drP, txP, prm, debug, mode);
            selected = 0;
        end
    end

    %% Drawing
%     Screen('FillRect', dpP.window, drP.grey);

    % Left button
    if selected == 1
        Screen('FillRect', dpP.window, drP.midGrey, leftRect);
    else
        Screen('FillRect', dpP.window, drP.grey, leftRect);
    end
    Screen('FrameRect', dpP.window,  drP.black, leftRect);
    DrawFormattedText(dpP.window, 'Demo Mouse', 'center', 'center',  drP.white, [], [], [], [], [], leftRect);

    % Right button
    if selected == 2
        Screen('FillRect', dpP.window, drP.midGrey, rightRect);
    else
        Screen('FillRect', dpP.window, drP.grey, rightRect);
    end
    Screen('FrameRect', dpP.window,  drP.black, rightRect);
    DrawFormattedText(dpP.window, 'Treino', 'center', 'center',  drP.white, [], [], [], [], [], rightRect);

    % Skip arrow
    Screen('DrawText', dpP.window, '>> Avançar', skipRect(1), skipRect(2),  drP.white);

    % Desenha o cursor antes de atualizar a tela
    if cursorVisible
        Screen('FillOval', dpP.window, prm.cursorColor, ...
            [mx-prm.cursorRadius_px my-prm.cursorRadius_px mx+prm.cursorRadius_px my+prm.cursorRadius_px]);
    end

    Screen('Flip', dpP.window);
end
end

function [nTs, targetOri, modTimes, nStimsToReport, orderToReportSets] = getForagingDistributions(nTargets, nStims, nTrials, nBlocks, params)

% (a) Distribuição da quantidade de alvos: triangular (discreta) ao 
    %     redor da média nTargets, se nTargets \in (1, nStims)
        if nTargets > 1 && nTargets < nStims
            numTargetsPMF = [.25 .50 .25];
            targetsSpace  = nTargets + [-1 0 1];
            nTs = randsample(targetsSpace, nBlocks*nTrials, true, numTargetsPMF);
            nTs = reshape(nTs, [nBlocks nTrials]);
        else
            nTs = nTargets*ones(nBlocks, nTrials);
        end

    % (b) Distribuição da orientação-alvo: uma por bloco, uniforme nas
    %     orientações de allOri
        targetOri = randsample(params.allOri, nBlocks, true);

    % (c) Distribuição do instante de modificação: após quantos alvos fixados 
    %     será apresentado o ruído rosa com orientação. A distribuição é
    %     tal que, em média, 3/4 do total de estímulos tenham sido vistos
        modeDistr = floor(3*nStims/4);
        targetModTimePMF = robust_beta_pmf(nStims, modeDistr, 'peakness', 1.5);
        auxModTimes = randsample(1:nStims, nBlocks*nTrials, true, targetModTimePMF);
        auxModTimes = reshape(auxModTimes, nBlocks, nTrials);
        modTimes = auxModTimes - 1;

    % (d) Distribuição da quantidade de estímulos cuja orientação deve ser
    %     reportada por trial: uniforme entre minToReport e maxToReport,
    %     e em propTrialsPSA*nTrials trials será verificado o efeito de
    %     atenção pré-sacádica
        nStimsToReport = round((params.maxToReport - params.minToReport)*rand(nBlocks,nTrials) + params.minToReport);
        reportPostModFix = rand(nBlocks, nTrials) <= params.propTrialsPSA;
        % Guarda a quantidade de estímulos a serem reportados que não são 
        % alvos da fixação 
        nStimsToReport = nStimsToReport - reportPostModFix; 
        % Garante que, se houver 2 ou mais estímulos (não alvos de fixação) 
        % a serem reportados, haja pelo menos um já visto e um não visto
        nStimsBase = nStimsToReport >= 2;
        % Tira 2 dos maiores (ou iguais a) 2, para distribui-los, e 0 dos
        % menores (i.e., de tamanho 1). Do que sobra, joga aleatoriamente
        % para cada lado
        aux = nStimsToReport-2*nStimsBase;
        nStimsAddPre  = round(rand(size(aux)).*aux);
        nStimsAddPost = (nStimsToReport - 2*nStimsBase) - nStimsAddPre;
        reportStimsFromSets = cat(1, permute(nStimsAddPre + nStimsBase,[3,2,1]), permute(nStimsAddPost + nStimsBase,[3,2,1]));
        nStimsToReport = cat(1, reportStimsFromSets(1,:,:), permute(reportPostModFix, [3, 2, 1]), reportStimsFromSets(2,:,:));
        % Em cada trial, o sujeito deve reportar dos estímulos vistos
        % anteriormente, do atual ou dos não vistos em ordem aleatória
        [~, orderToReportSets] = sort(rand(3, nTrials, nBlocks), 1);
end

function runForaging(tkP, dpP, drP, txP, params, debug, mode)

    modeMap = containers.Map({'cursor', 'training', 'experiment'}, 1:3);
    mode = modeMap(mode);

    
    leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4};



%     cursorVisible = false;
    if mode == 1
        ShowCursor(dpP.window);
%         cursorVisible = true;
        params.minFixTime2 = params.minCursorFixTime2;
        params.minFixTime3 = params.minCursorFixTime3;
        tkP.nBlocks = 1;
        tkP.nTrials = 5;
        % No treino, todas as orientações são alvo uma vez
    elseif mode == 2
        L = numel(params.allOri);
        tkP.nBlocks = L;
        tkP.nTrials = 5;
    end

    %% 2) Define as distribuições das condições dos trials e dos blocos
        [nTs, targetOri, modTimes, nStimsToReport, orderToReportSets] = getForagingDistributions(tkP.nTargets, tkP.nStims, tkP.nTrials, tkP.nBlocks, params);

    if mode == 2,  targetOri = randperm(L); end
    % (d) Matrizes com orientação e centros dos estímulos. Sem os alvos, há uma 
    %     quantidade aleatória de estímulos diagonais e cardinais
        stimCenters = zeros(2, tkP.nStims, tkP.nTrials, tkP.nBlocks);
        orientation = zeros(tkP.nStims, tkP.nTrials, tkP.nBlocks);
        for b=1:tkP.nBlocks
            auxAllOri = params.allOri;
            auxAllOri(auxAllOri == targetOri(b)) = [];
            orientation(:,:,b) = params.allOri(randi(length(auxAllOri), tkP.nStims, tkP.nTrials));
        end

        
%% 4) Qualidades das cruzes de fixação
        crossSize_px = dva2pix(params.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, params.crossSize_dva);
        fixCenters  = zeros(2, tkP.nTrials, tkP.nBlocks);
    
%% 5) Qualidades dos ruídos
    % (c) Define a matriz que contém o ruído e o centro dos retângulos a serem
    %     plotados
        noiseMatrix = zeros(txP.gabor.size_px, tkP.nStims*txP.gabor.size_px);
        oriPinkMatrix = zeros(txP.gabor.size_px, tkP.nStims*txP.gabor.size_px);
    
        noiseCenters = [0:txP.gabor.size_px:(tkP.nStims-1)*txP.gabor.size_px; zeros(1, tkP.nStims)] + txP.gabor.size_px/2;
    
    % (d) Centros de cada patch de ruído
        baseRect = [0 0 txP.gabor.size_px txP.gabor.size_px];
        srcRects = CenterRectOnPointd(repmat(baseRect, [tkP.nStims,1])', noiseCenters(1,:), noiseCenters(2,:));

        minDist_px   = dva2pix(params.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, params.minDist_dva);
        minFixDist1 = txP.gabor.size_px*params.fixDistFactor1;
    

    
%% 7) Localização dos estímulos (pré-computada para todos trials)
        for b=1:tkP.nBlocks
        for i=1:tkP.nTrials
    % (a) Cria os centros dos retângulos que contêm os Gabores e a fixação
            [currFixCenter, currStimCenter] = getStimLocations(dpP.ellipseProps, tkP.nStims, minDist_px, params.gridShape, params.randomize);
            fixCenters(:, i, b)     = currFixCenter;
            stimCenters(:, :, i, b) = currStimCenter;
    % (b) Atribui orientação não nula para apenas A(i) dos Gabores
            orientation(randperm(tkP.nStims, nTs(b, i)), i, b) = targetOri(b);
        end
        end
    % c) Distância mínima para considerar fixação em alvo pós-modificação
    minFixDist3 = txP.gabor.size_px*params.fixDistFactor3;
    try
        Screen('TextFont', dpP.window, 'Times');
        for b = 1:tkP.nBlocks
            restartBlock = false;
            % (-a) Mensagem com orientação a ser procurada

            Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
            Screen('DrawTexture', dpP.window, txP.exampleNoise.tex, [], [], targetOri(b), [], [], [], []);
            Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
            Screen('DrawTexture', dpP.window, txP.exampleGabor.tex, [], [], targetOri(b), [], [], [], [], [], txP.exampleGabor.props);
            Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTextures', dpP.window, txP.exampleBlob.tex, [], [], targetOri(b), [], [], [0 0 0 1]', [], [], txP.exampleBlob.props); 

            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('TextSize', dpP.window, 45);
            countBlockText  = sprintf('Bloco %d de %d', b, tkP.nBlocks);
            DrawFormattedText(dpP.window, countBlockText, 'center', dpP.screenRes.height*.1, drP.black);

            Screen('TextSize', dpP.window, 60);
            targetText_1    = 'Orientação dos alvos:';
            DrawFormattedText(dpP.window, targetText_1, 'center', dpP.screenRes.height*.2, drP.black);

            targetText_2    = params.allOriName{params.allOriMap(targetOri(b))};
            Screen('TextStyle', dpP.window, 1); Screen('TextSize', dpP.window, 80);
            DrawFormattedText(dpP.window, targetText_2, 'center', dpP.screenRes.height*.3, drP.black);
                                                            

            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('TextStyle', dpP.window, 0);
            Screen('TextSize', dpP.window, 60);
            proceedText = 'Pressione ESPAÇO para prosseguir';
            DrawFormattedText(dpP.window, proceedText, 'center', dpP.screenRes.height*.75, drP.black);

            Screen('TextSize', dpP.window, 24);

            Screen('Flip', dpP.window);

            while true
                [~,~,keyCode] = KbCheck;
                if keyCode(spaceKey), break; end
                WaitSecs(.01);
            end
            blockOnset = Screen('Flip', dpP.window);
            if debug == 0 && mode >= 2
                if mode == 2, prefix = 'trn'; else, prefix = 'exp'; end
                Eyelink('Message', sprintf('%s_Block_onset_%d_%d', prefix, b, blockOnset));
            end

        for i = 1:tkP.nTrials
        % (a) Cria os retângulos de destino
            dstRects = CenterRectOnPointd(repmat(baseRect, [tkP.nStims,1])', stimCenters(1,:, i, b), stimCenters(2,:, i, b));
        
        % i. Deixa de registrar até StartRecording
            if debug == 0 && mode >= 2
                Eyelink('SetOfflineMode');
            end
            WaitSecs(.1);

        
            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            % ii. Faz drift correction
            if debug == 0 && mode >= 2
                EyelinkDoDriftCorrection(el);
                WaitSecs(.1);
        % iii. Desenha, na tela do Host PC, linhas indicando a orientação de cada estímulo

                Eyelink('Command', 'clear_screen 0');
                for j=1:size(dstRects, 2)
                    hostRect = dstRects(:,j);
                    cx = mean(hostRect([1 3]));
                    cy = mean(hostRect([2 4]));
                    hostRect = round(hostRect);
                    radius = (hostRect(3)-hostRect(1))/2;
                    theta = orientation(j, i, b); lineLen = radius * 0.8;
                    dx = lineLen * sind(theta); dy = lineLen * cosd(theta);
                    Eyelink('Command', 'draw_line %d %d %d %d 12', round(cx-dx), round(cy-dy), round(cx+dx), round(cy+dy));
                end
            end
        
        % (a) Desenha a cruz de fixação
            xFix = [-crossSize_px/2 crossSize_px/2 0 0];
            yFix = [0 0 -crossSize_px/2 crossSize_px/2];
            fixCoords = [xFix; yFix];
        
            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawLines', dpP.window, fixCoords, params.lineWidth_px, drP.white, fixCenters(:, i, b)', 2);
        
        % (b) Cria as texturas de ruído, sem desenhá-las
            auxNoiseMatrix = butterFilter(pinkNoise(txP.gabor.size_px, tkP.nStims*txP.gabor.size_px), txP.noiseLoCutFreq, txP.noiseHiCutFreq);
            noiseMatrix(:,:) = params.noiseAlpha*(params.noiseContrast*rescale(auxNoiseMatrix, -params.noiseAmplitude, params.noiseAmplitude))+drP.grey;
            noiseTex = Screen('MakeTexture', dpP.window, noiseMatrix, [], [], [], [], []);

        % (c) Cria as texturas dos ruídos com orientação, sem desenhá-las
            for j=1:tkP.nStims
                colRange = ((j-1)*txP.gabor.size_px+1):(j*txP.gabor.size_px);
                oriPinkMatrix(:,colRange) = ApplyOriFilter(txP.oriFilter', txP.OFsize, auxNoiseMatrix(:,colRange));
            end

            oriPinkTex = Screen('MakeTexture', dpP.window, oriPinkMatrix);

            FPonset = Screen('Flip', dpP.window);
        
        % iv. Inicia o registro da sessão
        if debug == 0 && mode >= 2
            Eyelink('StartRecording');
            Eyelink('Command', 'record_status_message "TRIAL %d/%d"', i, tkP.nTrials);
            Eyelink('Message', sprintf('FP_onset_%d', FPonset));
        end
        
        % v. Não avança de tela até que os olhos estejam na cruz de fixação
            fixCenter = fixCenters(:,i, b);
            if debug > 0 && mode > 1
                KbPressWait;
            else
                while true
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
                    end
                    if check
                        % Se os olhos estiverem perto da cruz por tempo
                        % suficiente, prossegue
                        if vecnorm([x_gaze; y_gaze] - fixCenter) <= minFixDist1
                            if (GetSecs - FPonset) >= params.minFixTime1, break; end
                        % Se estiverem distantes, reinicia a contagem
                        elseif vecnorm([x_gaze; y_gaze] - fixCenter) > minFixDist1
                            FPonset = GetSecs;
                        end
                    end
                    WaitSecs(.01);
                    [~, ~, keyCode] = KbCheck;
                    if keyCode(escapeKey)
                        restartBlock = true;
                    end
                end
            end
            
            vbl = Screen('Flip', dpP.window);
            
        % (d) Desenha os ruídos, puramente opacos
            Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
            Screen('DrawTextures', dpP.window, noiseTex, srcRects, dstRects, orientation(:, i, b), [], [], [], []);
        
        % (e) Desenha os gratings, somando ambos os sinais
            Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
            Screen('DrawTextures', dpP.window, txP.gabor.tex, [], dstRects, orientation(:, i, b), [], [], [], [], [], txP.gabor.props);
        
        % (f) Desenha a abertura gaussiana
            Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects, orientation(:, i, b), [], [], [0 0 0 1]', [], [], txP.blob.props);                                                             
            Screen('Close', noiseTex);
        
            vbl = Screen('Flip', dpP.window, vbl + 0.5 * dpP.ifi);

        % vi. Registra o momento de início do trial
            trialOnset = vbl;
            if debug == 0  && mode >= 2
                Eyelink('Message', sprintf('trial_onset_%1d', trialOnset));
                Eyelink('Command', 'record_status_message "TRIAL %d', i);
            end

        % vii. O sujeito deve visitar exatamente Ns(i) estímulos
        %      antes de ocorrerem as modificações pré-sacádicas
            counter = 0;                % Número de estímulos visitados
             flag = zeros(1, tkP.nStims); % Quantas vezes cada estímulo foi visitado

            if debug > 0 && mode > 1
                seenIdx = sort(randsample(tkP.nStims,modTimes(b, i)));
                flag(seenIdx) = 1;
                currIdx = randsample(seenIdx, 1);
                KbPressWait;
            else
                while counter < modTimes(b, i)
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
                    end

                    if check
                        % Se os olhos estiverem suficientemente próximos de
                        %  um alvo não antes visto, incrementa o contador
                        isCurrTarget = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, i, b)) <= minFixDist1;
                        if any(isCurrTarget)
                            if (GetSecs - trialOnset) >= params.minFixTime2
                                if any(isCurrTarget & (flag == 0))
                                    counter = counter + 1;
                                end
                                flag = flag + isCurrTarget;
                            end
                        % Se estiverem distantes, reinicia o tempo
                        else
                            % disp('Sem fixação')
                            trialOnset = GetSecs;
                        end
                    end
                    WaitSecs(0.001)
                end
                currIdx = find(isCurrTarget);
                seenIdx = find(flag ~= 0);
            end
            notSeenIdx = find(flag == 0);


        % (g) Desenha ruído orientado em todos os estímulos menos o atual
            % (as linhas comentadas servem para mudar apenas os não vistos
            blinkIdx = setdiff(1:tkP.nStims, currIdx);
            Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
%             Screen('DrawTextures', window, oriPinkTex, srcRects(:,notSeenIdx), dstRects(:,notSeenIdx), orientation(notSeenIdx, i, b), [], [], [], []);
            Screen('DrawTextures', dpP.window, oriPinkTex, srcRects(:,blinkIdx), dstRects(:,blinkIdx), orientation(blinkIdx, i, b), [], [], [], []);
            Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
%             Screen('DrawTextures', window, blobTex, [], dstRects(:,notSeenIdx), orientation(notSeenIdx, i, b), [], [], [0 0 0 1]', [], [], blobProps);
            Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects(:,blinkIdx), orientation(blinkIdx, i, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
            Screen('Close', oriPinkTex);
    
        % viii. Registra os tempos de início e fim da 3a fase
            updateStimOnset = Screen('Flip', dpP.window);
            if debug == 0 && mode >= 2
                Eyelink('Message', sprintf('update_onset_%1d', updateStimOnset));
                updateStimOffset = Screen('Flip', dpP.window, updateStimOnset + params.pinkNoiseDur);
                Eyelink('Message', sprintf('update_offset_%1d', updateStimOffset));
            else
                updateStimOffset = Screen('Flip', dpP.window, updateStimOnset + params.pinkNoiseDur);
            end

        % ix. Verifica se há alguma fixação de duração mínima em estímulo 
        %     numa janela pós-modificação
            fixOnset = updateStimOffset; currIdx = [];

            if debug == 0 || mode == 1
                while GetSecs - updateStimOffset < params.postModDur
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
                    end

                    if check
                        isCurrTarget = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, i, b)) <= minFixDist3;
                        if any(isCurrTarget)
                            if (GetSecs - fixOnset) >= params.minFixTime3
                                currIdx = find(isCurrTarget);
                                disp(['Trial ' num2str(i) ': Visitou o alvo ' num2str(currIdx) ' pós-modificação'])
                                break;
                            end
                        else
                            fixOnset = GetSecs;
                        end
                    end
                    WaitSecs(.001)
                end
            end

            if ~isempty(currIdx) && ismember(currIdx, notSeenIdx)
                notSeenIdx(notSeenIdx == currIdx) = [];
            end
        
        % x. Desenha placeholders para os estímulos aleatoriamente, obedecendo
        %    a identificação deles como visitados, atual e não visitados
            allTargets = nan(1,tkP.nStims); allColors2 = drP.allColors;
            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            seenAux = datasample(setdiff(seenIdx, currIdx), min(length(seenIdx), nStimsToReport(1, i, b)), 'Replace', false);
            currAux = []; if nStimsToReport(2, i, b) == 1, currAux = currIdx; end
            notSeenAux = datasample(notSeenIdx, min(length(notSeenIdx), nStimsToReport(3, i, b)), 'Replace', false);
            orderToReportStimsCell = {seenAux, currAux, notSeenAux};
            orderToReportStims = [orderToReportStimsCell{orderToReportSets(1, i, b)} orderToReportStimsCell{orderToReportSets(2, i, b)} orderToReportStimsCell{orderToReportSets(3, i, b)}];
            
            if debug == 0  && mode >= 2
                Eyelink('Message', sprintf('answer_onset_%1d', GetSecs));
            end
            feedback = zeros(1, length(orderToReportStims));
            for j=1:length(orderToReportStims)
                KbReleaseWait;
                rectColors = allColors2; rectColors(:, orderToReportStims(j)) = drP.orange;
            
                rectPW = drP.allPW; rectPW(:, orderToReportStims(j)) = params.pW2;
                currTarget = false;
                abort  = false;
            
                while ~abort
                    [keyIsDown, ~, keyCode] = KbCheck;
                    if keyIsDown
                        if keyCode(leftKey)
                            currTarget = false;
                        elseif keyCode(rightKey)
                            currTarget = true;
                        elseif keyCode(spaceKey)
                            abort = true;
                            allColors2(:, orderToReportStims(j)) = drP.midGrey;
                        end
                    end
            
                    allTargets(orderToReportStims(j)) = currTarget;
            
                    answerOffset = foragingFlip(dpP.window, stimCenters(:, :, i, b), dstRects, orderToReportStims, txP.gabor.size_px, rectColors, allTargets, targetOri(b), rectPW);
            
                end
            end
            
            if debug == 0 && mode >= 2
                Eyelink('Message', sprintf('answer_offset_%1d', answerOffset));
            end
            if mode < 3
                feedback = (orientation(:,i,b) == targetOri(b))' + allTargets;
                feedback(rem(feedback,2) == 0) = 2; feedback(rem(feedback,2) == 1) = 0;
                feedback = feedback/2;
                foragingFlip(dpP.window, stimCenters(:, :, i, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW, feedback, drP.red, drP.green);
            else
                foragingFlip(dpP.window, stimCenters(:, :, i, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW);
             end
            WaitSecs(1);

        % xi. Interrompe o registro, pois ou a tela será atualizada ou
        %     acabaram os trials
            
            if debug == 0 && mode >= 2
                Eyelink('Message', sprintf('trial_offset_%1d', GetSecs));
                Screen('Flip', dpP.window);
                WaitSecs(0.1);
                Eyelink('SetOfflineMode');
                Eyelink('StopRecording');
            else 
                Screen('Flip', dpP.window);
                WaitSecs(0.5);
            end
        end
        end
    catch
        psychrethrow(psychlasterror);
        cleanup;
    end

    if mode == 1, HideCursor(dpP.window); end
end

function gabor = foragingGabor(window, screenRes, nStims, params, monitorW_mm)
    % (a) Tamanho dos gabores, distância entre eles, em px, e cor do fundo
        gaborSize_px = dva2pix(params.screenDist, monitorW_mm/10, screenRes.width, params.gaborSize_dva);
        gaborBackgroundOffset = [0 0 0 0];

    % (b) Cria a textura em si
        gaborTex = CreateProceduralSineGrating(window, gaborSize_px, gaborSize_px, gaborBackgroundOffset,[], 0.5);
    % (c) Propriedades do Gabor
        gaborPPD = gaborSize_px/params.gaborSize_dva;
        gaborFreq = params.gaborFreq_cpd/gaborPPD;    % Frequência, ciclos por pixel (cpp)
        gaborSigma = gaborSize_px / 6;
        
        gaborProps0 = [params.gaborPhase, gaborFreq, params.gaborAlpha*params.gaborAmplitude*params.gaborContrast, 0];
        gaborProps = repmat(gaborProps0', 1, nStims);

        gabor.size_px = gaborSize_px;
        gabor.ppd     = gaborPPD;
        gabor.props   = gaborProps;
        gabor.sigma   = gaborSigma;
        gabor.tex     = gaborTex;
end

function noise = foragingNoise(window, noiseSize, nStims, grey, params, loCut, hiCut)
        auxNoiseMatrix = pinkNoise(noiseSize, nStims*noiseSize);
        if nargin >= 6
            auxNoiseMatrix = butterFilter(auxNoiseMatrix, loCut, hiCut);
        end
        noiseMatrix = params.noiseAlpha*(params.noiseContrast*rescale(auxNoiseMatrix, -params.noiseAmplitude, params.noiseAmplitude))+grey;
        noiseTex = Screen('MakeTexture', window, noiseMatrix, [], [], [], [], []);
        noise.size_px = [noiseSize nStims*noiseSize];
        noise.tex     = noiseTex;        
end

function blob  = foragingBlob(window, blobSize, grey, params)
        blobColorOffset = [grey grey grey 0];
        [blobTex, ~] = CreateProceduralGaussBlob(window, blobSize, blobSize, blobColorOffset, 1, 1);
        blobSigma = blobSize/6;
        blobProps = [params.blobContrast, blobSigma, params.blobAspect, 0]';

        blob.size_px = blobSize;
        blob.props   = blobProps;
        blob.sigma   = blobSigma;
        blob.tex     = blobTex;
end
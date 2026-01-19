function ForagingGabors(nTrials, nStims, nBlocks, nTargets, options)
% A tarefa consiste em nTrials trials, cada um com nStims estímulos e 
% quantidade média de nTargets alvos, espalhados pseudoaleatoriamente na 
% tela. A tarefa do sujeito é encontrar, dentre gabores de alta frequência 
% com orientações variadas, alvos com uma orientação específica, apenas com
% o movimento ocular. Em determinado momento, durante essa busca ocular, os
% estímulos ainda não visitados são substituídos por um ruído rosa de mesma 
% orientação durante o período pré-sacádico. Por fim, o sujeito deve
% reportar se viu ou não alvo nas posições indicadas

%% Parâmetros que variam conforme o participante:
% (a) Olho dominante, (b) mediana do tempo de fixação, (c) parâmetros do ruído rosa;
% (d) tecla de seleção de alvo, (e) ordem da orientação-alvo.

%% Convenção:
% Os objetos estão indexados na ordem (nStims, nTrials, nBlocks) quando têm
% 3 dimensões ou como (nBlocks, nTrials) quando têm 2 dimensões (i.e.,
% nTrials é sempre a 2a dimensão)
    
    arguments
        nTrials {mustBeNumeric}   = 20;
        nStims {mustBeNumeric}    = 8;
        nBlocks {mustBeNumeric}   = 15;
        nTargets {mustBeNumeric}  = ceil(nStims/2);
        options.mode string       = 'experiment' % 'experiment', 'debug' ou 'debugTV'
    end
    
    rng('shuffle');
    cleanup

    % Verifica se o programa será executado no modo debug ou não
    debug = 0;
    if strcmp(options.mode, 'debug'), debug = 1; elseif strcmp(options.mode, 'debugTV'), debug = 2; end
    clear options;

    if debug > 0
        nTrials  = 2;
        nStims   = 8;
        nBlocks  = 2;
        nTargets = ceil(nStims/2);
    end

%% 1) Importa os parâmetros para o experimento
    params = foragingParams;
    params.filePath = fileparts(mfilename('fullpath'));
    params.edfExtension = '.edf';
    
%% 2) Inicializa PTB e EyelinkToolBox
    % i. Caixa de diálogo
    edfFile = ''; targetKey = 'left';
    if debug == 0
        repeat = true;
        while repeat
            repeat = false;
            prompt = {'\fontsize{14} Sujeito', '\fontsize{14} Número da sessão', '\fontsize{14} Olho dominante (E ou D)', '\fontsize{14} Seta para alvo (E ou D)'};
            dlg_title = 'Info. sessão';
            def = {'00', '01', 'D', 'D'};
            options.Resize = 'on'; options.Interpreter = 'tex';
            answer = inputdlg(prompt, dlg_title, [1 40; 1 40; 1 40; 1 40], def, options);
    
            % Se a caixa de diálogo for fechada, encerra a sessão
            if isempty(answer), fprintf('Sessão cancelada pelo usuário\n'); cleanup; return; end
    

            % Recomeça se nome muito longo
            edfFile = [answer{1} '_' answer{2}];
            if length(edfFile) > 8, fprintf('ERRO: O nome do arquivo excede 8 caracteres!\nTente novamente\n'); cleanup; repeat = true; end
    
            % Recomeça se campo de olho dominante mal preenchido
            domEye = upper(answer{3});
            if ~ismember(domEye, {'E', 'L', 'D', 'R'}), fprintf('ERRO: Preencha o campo de olho dominante com E para esquerdo ou D para direito.\nTente novamente\n'); cleanup; repeat = true; end
            targetKey = upper(answer{4});
            if ~ismember(targetKey, {'E', 'L', 'D', 'R'})
                fprintf('ERRO: Preencha o campo de deta para alvo com E para esquerda ou D para direita.\nTente novamente\n'); 
                cleanup; 
                repeat = true;
            else
                if ismember(targetKey, {'E', 'L'})
                    targetKey = 'left';
                elseif ismember(targetKey,{'D', 'R'})
                    targetKey = 'right';
                end
            end
        end
    end
    clear repeat prompt dlg_title def options answer
    
    % (a) Configurações antes de abrir a tela
        
    FlushEvents;
    PsychDefaultSetup(2);   % Dispensa KbName('UnifyKeyNames')
    Screen('Preference','TextAlphaBlending',1);
    Screen('Preference', 'SuppressAllWarnings', 1);

    if debug == 0
        Screen('Preference', 'SyncTestSettings', [], [], 0.2, 10);
    else
        Screen('Preference', 'SkipSyncTests', 2);
        Screen('Preference','VisualDebugLevel', 0);
        Screen('Preference', 'Verbosity', 0);
    end
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','General','UseFastOffscreenWindows');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

    leftKey   = KbName('LeftArrow'); rightKey  = KbName('RightArrow');
    spaceKey  = KbName('space');     escapeKey = KbName('ESCAPE');
    rKey      = KbName('r');
    keys = {leftKey, rightKey, spaceKey, escapeKey, rKey};

    % (b) Escolhe a tela em que haverá o desenho e define algumas cores
    screenNumber = max(Screen('Screens'));
    black = BlackIndex(screenNumber);
    white = WhiteIndex(screenNumber);
    grey = white / 2;
    if debug == 1, bgColor = black; else, bgColor = grey; end
    [window, winRect] = PsychImaging('OpenWindow', screenNumber, bgColor, [], 32, 2, [], [],  kPsychNeed32BPCFloat);
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % (c) Obtém propriedades da tela
    [monitorW_mm, ~] = Screen('DisplaySize', screenNumber); % Tamanho   da tela em mm
    screenRes = Screen('Resolution', screenNumber);         % Resolução da tela em px
    ifi = Screen('GetFlipInterval', window);

    % (d) Obtém o centro e os semieixos da elipse
    ellipseProps = [screenRes.width/2, screenRes.height/2];
    ellipseProps = [ellipseProps ellipseProps(1)*params.ellipseToScreenRatio(1) ellipseProps(2)*params.ellipseToScreenRatio(2)];

    % (e) Define outras cores de interesse e espessura das canetas para
    %     desenhar algumas das figuras
    whiteGrey = (white+grey)/2;
    blackGrey = (black+grey)/2;
    orange   = white*params.orange;
    blue     = white*params.blue;
    darkBlue = white*params.darkBlue;
    red     = white*params.red;
    green   = white*params.green;
    params  = rmfield(params, {'blue', 'darkBlue', 'green', 'orange', 'red'});
    
    allColors = white*ones(3, nStims);
    allPW     = params.pW1*ones(1,nStims);

    % i. Inicializa Eyelink
    dummymode = 0; if debug > 0, dummymode = 1; disp('Debug: O EyeLink não será inicializado'); end
    EyelinkInit(dummymode);

    % ii. Abre o arquivo .edf
    el = {};
    if debug == 0
        failOpen = Eyelink('OpenFile', edfFile);
        % Interrompe se arquivo não abrir
        if failOpen ~= 0, fprintf('ERRO: Não foi possível criar o arquivo %s\n', edfFile); cleanup; return; end
    
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
        ListenChar(-1); % Para o teclado não escrever na command window

        clear failOpen domEye preambleText topPriorityLevel
    else
        Eye = [];
    end

    clear dummymode leftKey rightKey spaceKey escapeKey rKey screenNumber

%% 3) Obtém as texturas de Gabor, abertura gaussiana e do ruído
    if debug == 2, params.screenDist = 200; end

    % (a) Estímulos para a tarefa
    gabor = foragingGabor(window, screenRes, nStims, params, monitorW_mm);
    blob  = foragingBlob(window, gabor.size_px, grey, params);

    noiseHiCutFreq = params.noiseHiCutFreq_cpd/gabor.ppd; % Parâmetros do filtro para ruído
    noiseLoCutFreq = params.noiseLoCutFreq_cpd/gabor.ppd;

    [oriFilter, OFsize] = MakeOriFilter(gabor.size_px, params.aSigma, params.rSigma2); % Filtro de orientação

    % (b) Estímulos para a tela entre blocos
    exampleParams = params; exampleGaborFactor = 3;
    exampleParams.gaborFreq_cpd = params.gaborFreq_cpd/exampleGaborFactor;
    exampleParams.gaborSize_dva = params.gaborSize_dva*exampleGaborFactor;
    exampleGabor = foragingGabor(window, screenRes, 1, exampleParams, monitorW_mm);
    exampleBlob = foragingBlob(window, exampleGabor.size_px, grey, params);
    exampleNoise = foragingNoise(window, exampleGabor.size_px, 1, grey, params);

%% 4) Reúne todas as variáveis definidas até aqui para passá-las às demais funções
    fixQueue = params.minFixTime2*ones(1, max(params.fixTimeQueueSize, round(2.5*nStims)));

    displayProps.window       = window;
    displayProps.winCol       = grey;
    displayProps.winRect      = winRect;
    displayProps.winCenter    = [winRect(3)/2 winRect(4)/2];
    displayProps.monitorW_mm  = monitorW_mm;
    displayProps.screenRes    = screenRes;
    displayProps.ifi          = ifi;
    displayProps.ellipseProps = ellipseProps;
    clear exampleParams exampleGaborFactor bgColor window winRect monitorW_mm screenRes ifi ellipseProps

    drawProps.white     = white;
    drawProps.grey      = grey;
    drawProps.whiteGrey = whiteGrey;
    drawProps.blackGrey = blackGrey;
    drawProps.black     = black;
    drawProps.orange    = orange;
    drawProps.blue      = blue;
    drawProps.darkBlue  = darkBlue;
    drawProps.red       = red;
    drawProps.green     = green;
    drawProps.allColors = allColors;
    drawProps.allPW     = allPW;
    clear white whiteGrey blackGrey grey black orange blue darkBlue red green allColors allPW

    texProps.gabor          = gabor;
    texProps.blob           = blob;
    texProps.exampleGabor   = exampleGabor;
    texProps.exampleBlob    = exampleBlob;
    texProps.exampleNoise   = exampleNoise;
    texProps.noiseHiCutFreq = noiseHiCutFreq;
    texProps.noiseLoCutFreq = noiseLoCutFreq;
    texProps.oriFilter      = oriFilter;
    texProps.OFsize         = OFsize;
    clear gabor blob exampleGabor exampleBlob exampleNoise noiseHiCutFreq noiseLoCutFreq oriFilter OFsize

    taskProps.Eye       = Eye;
    taskProps.targetKey = targetKey;
    taskProps.el        = el;
    taskProps.edfFile   = edfFile;
    taskProps.fixQueue  = fixQueue;
    taskProps.nTargets  = nTargets;
    taskProps.nStims    = nStims;
    taskProps.nTrials   = nTrials;
    taskProps.nBlocks   = nBlocks;
    taskProps.keys      = keys;
    clear ans Eye targetKey el edfFile fixQueue nTargets nStims nTrials nBlocks keys

%% 5) Chama as funções para executar a tarefa
%     if debug ~= 1
    fakeLoadingScreen(taskProps, displayProps, drawProps, params);
%     else
%         Screen('FillRect', displayProps.window, drawProps.grey);
%     end
    [~, taskState] = menuScreen(taskProps, displayProps, drawProps, texProps, debug, params);

    warningStart = 'AVISO: Sessão ';
    myWarning = {};
    warnings = {'de treino não concluída.', 'encerrada sem treino ou experimento.', 'de treino concluída e salva.'; ...
                'experimental não concluída nem salva.', 'experimental não concluída, mas salva.', 'experimental concluída e salva com êxito.'};

    if sum(taskState(:)) == 0
        myWarning{end+1} = [warningStart warnings{1, 2}];
        disp(myWarning{end});
    else
        for i=1:size(taskState,1)
            if sum(taskState(i,:)) ~= 0
                colIdx = find(taskState(i,:) > 0, 1, "last");
                myWarning{end+1} = [warningStart warnings{i, colIdx}]; %#ok<AGROW> 
                disp(myWarning{end});
            end
        end
    end
    cleanup(displayProps.window);
end


function [tkP, taskState] = menuScreen(tkP, dpP, drP, txP, debug, prm, menuMode)
    if nargin < 7, menuMode = 'normal'; end
%% Interpretação da variável taskState:
    % Linhas: treino e experimento; 
    % Colunas: começou, salvou e concluiu
    taskState = [0 0 0;
                 0 0 0];

    Screen('TextFont', dpP.window, prm.textFont);
    Screen('TextSize', dpP.window, prm.textSizeNormal);
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Reposiciona o mouse
    cX = dpP.winCenter(1); cY = dpP.winCenter(2);
    SetMouse(cX, cY, dpP.window);
    if debug == 0,  HideCursor(dpP.window); end

    titleText = 'MENU'; titleMargin = cX*prm.titleMarginFactor;

    % Cria texturas de cursor e treino
    if strcmp(menuMode, 'normal')
        options     = {'cursor', 'training'};
        optionsName = {'Demo Mouse', 'Treino'};
        optionsMsg = {['Demonstração com o cursor (mouse), para\n'...
               'familiarização com a estrutura da tarefa.'],...
               ['Treino para o experimento, já com movimentos\n'...
               'oculares, e feedback ao fim de cada tentativa.\n\n'...
               'Recomendado antes de toda sessão experimental.'],...
               'Inicia a sessão experimental propriamente.'};
    else
        options     = {'staircase', 'cursor', 'training'};
        optionsName = {'Staircase', 'Demo Mouse', 'Treino'};
        optionsMsg = {['Ajusta os parâmetros do ruído conforme o\n'...
               'participante usando o método de staircase.'],...
               ['Demonstração com o cursor (mouse), para\n'...
               'familiarização com a estrutura da tarefa.'],...
               ['Treino para o experimento, já com movimentos\n'...
               'oculares, e feedback ao fim de cada tentativa.\n\n'...
               'Recomendado antes de toda sessão experimental.'],...
               'Inicia a sessão experimental propriamente.'};
    end
    parentDir = prm.filePath;
    iconsDir = cellfun(@(x) fullfile(parentDir, prm.imgFolder, [x prm.imgExtension]), options, 'UniformOutput', false);
    clear parentDir

    iconsTex = getMenuTex(dpP.window, iconsDir, drP.black);
    
    leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4};
    
    L = numel(options);

    %% Retângulos dos botões
    gap  = 80;
    btnW = min(300, (dpP.winRect(3)-(L+1)*gap)/L); btnH = 300;
    totalW = L*btnW + (L-1)*gap;
    leftX0 = cX - totalW/2;
    leftX  = leftX0 + (0:L-1) * (btnW + gap);

    btnRects = zeros(4, L);
    btnRects(1,:) = leftX; btnRects(2,:) = cY - btnH/2;
    btnRects(3,:) = leftX + btnW; btnRects(4,:) = cY + btnH/2;
    
    prm.upFrac = 2/3;
    RectSplit = btnRects(2,:)+btnH*prm.upFrac;
    upRects = btnRects; upRects(4,:) = RectSplit; downRects = btnRects; downRects(2,:) = RectSplit;
    uRH = RectHeight(upRects(:,1)'); uRW = RectWidth(upRects(:,1)');
    uRHW = min(uRH, uRW);
    for i=1:L
        [cXaux, cYaux] = RectCenter(upRects(:,i));
        upRects(:,i) = CenterRectOnPoint([0 0, uRHW, uRHW]*prm.iconScaleFactor, cXaux, cYaux);
    end
    clear cXaux cYaux

    %% Retângulo de info adicional
    detailRect = CenterRectOnPoint([0 0 2*prm.btnW (1.2/3)*prm.btnH], cX, 1.75*cY);
    drawBox = false;
    
    %% Seta para sair do menu/avançar
    skipText = 'Avançar';
    skipRect = Screen('TextBounds', dpP.window, skipText); %[dpP.winRect(3)-150, dpP.winRect(4)-60, 110, 40];
    skipRect(1) = dpP.winRect(3)-150; skipRect(2) = dpP.winRect(4)-60;

    arrow = [-25 -8;
              20 -8;
              25  0;
              20  8;
             -25  8;
             -20  0];
    
    scale = 3.5;
    arrow = arrow*scale + [skipRect(1) + skipRect(3)/2 skipRect(2) + skipRect(4)/2];
    clear lRectSplit rRectSplit scale
    
    %% Estados
    if (L+1)/2 == floor((L+1)/2)
        firstSelectL = (L+1)/2;
        firstSelectR = (L+1)/2;
    else
        firstSelectL = floor((L+1)/2);
        firstSelectR = ceil((L+1)/2);
    end
    selected = 0;
    [prevMx, prevMy, prevBut] = GetMouse(dpP.window);
    prevMouse = [prevMx, prevMy, prevBut];
    
    isMouseMostRecent = false;
    keyWasDown = false;
    lastActionTime = 0; lastRepeatTime = -Inf;
    doAction = false;
    mode = '';

    while true

        % Lê posição do mouse
        [mx, my, buttons] = GetMouse(dpP.window);
        mouseMoved = norm([mx my] - prevMouse(1:2)) > 1;
        prevMouse = [mx my buttons];
        cursorVisible = false;
    
        if mouseMoved || isMouseMostRecent
            cursorVisible = true;
            isMouseMostRecent = true;
            selected = 0;
            for i=1:L
                if IsInRect(mx, my, btnRects(:,i))
                    selected = i;
                    break;
                end
            end
            if inpolygon(mx, my, arrow(:,1), arrow(:,2))
                selected = L+1;
            end
        end
    
        %% Desenho

        Screen('TextSize', dpP.window, prm.textSizeTitle);
        Screen('DrawText', dpP.window, titleText, titleMargin, titleMargin,  drP.blackGrey);
        Screen('TextSize', dpP.window, prm.textSizeNormal);
        if selected == 0, drawBox = false; end

        % Desenha os retângulos (e seta) e textos e texturas
        Screen('FillRect', dpP.window, drP.blue, btnRects);
        Screen('FillPoly', dpP.window, drP.blue, arrow);

        Screen('DrawTextures', dpP.window, iconsTex, [], upRects);

        Screen('TextStyle', dpP.window, 1);
        for i=1:L
            DrawFormattedText(dpP.window, optionsName{i}, 'center', 'center',  drP.black, [], [], [], [], [], downRects(:,i)');
        end
        Screen('TextSize', dpP.window, prm.textSizeNormalish);
        Screen('DrawText', dpP.window, skipText, skipRect(1), skipRect(2),  drP.black);

        % Se houver algum retângulo selecionado, adiciona contorno
        if selected ~= 0
            drawBox = true;
            msg = optionsMsg{selected};
            if selected < L+1
                Screen('FrameRect', dpP.window,  drP.darkBlue, btnRects(:, selected), prm.pW2);
            else
                Screen('FramePoly', dpP.window, drP.darkBlue, arrow, prm.pW2);
            end
        end
        Screen('TextStyle', dpP.window, 0);

        if drawBox
            Screen('FillRect', dpP.window, drP.whiteGrey, detailRect);
            Screen('TextSize', dpP.window, prm.textSizeNormalish);
            DrawFormattedText(dpP.window, msg, 'center', 'center',  drP.blackGrey, [], [], [], [], [], detailRect);
            Screen('TextSize', dpP.window, prm.textSizeNormal);
        end

        % Desenha o cursor antes de atualizar a tela
        if cursorVisible
            Screen('FillOval', dpP.window, prm.cursorColor, ...
                [mx-prm.cursorRadius_px my-prm.cursorRadius_px mx+prm.cursorRadius_px my+prm.cursorRadius_px]);
        end
    
        Screen('Flip', dpP.window);
    
        % Lê inputs do teclado e clique do mouse
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown || any(buttons)
            if keyIsDown
                if ~keyWasDown
                    lastActionTime = GetSecs;
                    doAction = true;
                elseif GetSecs - lastActionTime > prm.repeatDelay
                    if GetSecs - lastRepeatTime > prm.repeatRate
                        doAction = true;
                        lastRepeatTime = GetSecs;
                    end
                end
            else
                doAction = true;
            end

            if doAction
                doAction = false;
                if keyCode(leftKey)
                    if selected == 0
                        selected = firstSelectL;
                    else
                        selected = max(1, selected - 1);
                    end
%                     cursorVisible = false;
                    isMouseMostRecent = false;
                elseif keyCode(rightKey)
                    if selected == 0
                        selected = firstSelectR;
                    else
                        selected = min(L+1, selected + 1);
                    end
%                     cursorVisible = false;
                    isMouseMostRecent = false;
                elseif keyCode(spaceKey) || any(buttons)
                    KbReleaseWait;
                    while any(buttons)
                        [~, ~, buttons] = GetMouse(dpP.window);
                        WaitSecs(0.001); 
                    end
                    if selected ~=0
                        if selected < L+1
                            mode = options{selected};
                        else
                            mode = 'experiment';
                        end
                    end
                    if strcmp(mode,'training')
                        taskState(1,1) = 1;
                    elseif strcmp(mode,'experiment')
                        taskState(2,1) = 1;
                    end
                    if selected ~= 0
                        for i=1:L, Screen('Close', iconsTex(i)); end
                        Screen('Flip', dpP.window);
                        WaitSecs(0.05);
                        % Apenas importa a fila de fixações se a tarefa não for de cursor
                        if strcmp(mode, 'cursor')
                            [~, savedNdone] = runForaging(tkP, dpP, drP, txP, prm, debug, mode);
                        elseif strcmp(mode, 'staircase')
                            disp('Vai para staircase...')
                        else
                            [tkP, savedNdone] = runForaging(tkP, dpP, drP, txP, prm, debug, mode);
                        end
                        if strcmp(mode, 'training'),      taskState(1,2:3) = savedNdone; end
                        if strcmp(mode, 'experiment'), taskState(2,2:3) = savedNdone; end
                        % Se o experimento tiver sido concluído, não volta
                        % pro menu
                        if taskState(2,3) == 1
                            return;
                        else
                            iconsTex = getMenuTex(dpP.window, iconsDir, drP.black);
                            SetMouse(cX, cY, dpP.window);
                            selected = 0;
                        end
                    end
                elseif keyCode(escapeKey)
                    KbReleaseWait;
                    decision = pauseScreen(tkP, dpP, drP, prm, 'menu');
                    if strcmp(decision, 'quit')
                        Screen('Flip', dpP.window); 
                        return;
                    end
                end
            end
        else
            lastRepeatTime = -Inf;
        end
        keyWasDown = keyIsDown;
    end
end


function [iconsTex] = getMenuTex(win, iconsDir, col)
    % Cria texturas de cursor e treino
    iconsImg = cellfun(@(x) alphaRead(x), iconsDir, 'UniformOutput', false);
    imgSizes = cellfun(@(x) size(x), iconsImg, 'UniformOutput', false);
    auxImg   = cellfun(@(x) col*ones(x(1), x(2), 4), imgSizes, 'UniformOutput', false);
    iconsImg = cellfun(@(x, a) addAlpha(x,a), auxImg, iconsImg, 'UniformOutput', false);
    iconsTex = cellfun(@(x) Screen('MakeTexture', win, x), iconsImg);
end


function [tkP, savedNdone, fixCenters, stimCenters, orientation] = runForaging(tkP, dpP, drP, txP, prm, debug, mode)
        % O modo default é experiment
        if nargin < 7, mode = 'experiment'; end

        Screen('Flip', dpP.window);

        if strcmp(mode, 'experiment')
            fakeLoadingScreen(tkP, dpP, drP, prm, 'start');
        end

        modeMap = containers.Map({'cursor', 'training', 'experiment'}, 1:3);
        mode = modeMap(mode);

        suffix = prm.msg.suffix{mode};
        
        leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3}; escapeKey = tkP.keys{4}; rKey = tkP.keys{5};
        if strcmp(tkP.targetKey, 'right')
            targetKey = rightKey; nonTargetKey = leftKey;
        else
            targetKey = leftKey; nonTargetKey = rightKey;
        end
    
    %% 1) Se não estiver no modo experimental, ajusta o número de trials e blocos
        if mode < 3
            nBlocks = tkP.nBlocks;
            nTrials = tkP.nTrials;
            if mode == 1
                prm.minFixTime2 = prm.minCursorFixTime2;
                prm.minFixTime3 = prm.minCursorFixTime3;
                prm.postModDur  = prm.postModDurCursor;
                tkP.nBlocks = 1;
                tkP.nTrials = prm.nTrialsTrain;
            % No treino, todas as orientações são alvo uma vez
            elseif mode == 2
                L = numel(prm.allOri);
                tkP.nBlocks = L;
                tkP.nTrials = prm.nTrialsTrain;
            end
        end

    %% 2) Define as distribuições das condições dos trials e dos blocos
        nTrialsBuffered = tkP.nTrials + prm.nBufferTrials;
        [nTs, targetOri, modTimes, nStimsToReport, orderToReportSets] = getForagingDistributions(tkP.nTargets, tkP.nStims, nTrialsBuffered, tkP.nBlocks, prm);
        if mode == 2,  targetOri = prm.allOri(randperm(tkP.nBlocks)); end

    %% 3) Define matrizes usadas para os estímulos
        % (a) Matrizes com orientação e centros dos estímulos. Sem os alvos,
        %     há quantidade aleatória de estímulos diagonais e cardinais 
        %     (cf. (e) para adição de alvos)
        stimCenters = zeros(2, tkP.nStims, nTrialsBuffered, tkP.nBlocks);
        orientation = zeros(tkP.nStims, nTrialsBuffered, tkP.nBlocks);
        for b=1:tkP.nBlocks
            auxAllOri = prm.allOri;
            auxAllOri(auxAllOri == targetOri(b)) = [];
            orientation(:,:,b) = prm.allOri(randi(length(auxAllOri), tkP.nStims, nTrialsBuffered));
        end

        % (b) Tamanho e matrizes para as cruzes de fixação 
        crossSize_px = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.crossSize_dva);
        fixCenters  = zeros(2, nTrialsBuffered, tkP.nBlocks);

        % (c) Matrizes com os ruído e o centro dos retângulos a serem
        %     plotados (srcRect)
        noiseMatrix = zeros(txP.gabor.size_px, tkP.nStims*txP.gabor.size_px);
        oriPinkMatrix = zeros(txP.gabor.size_px, tkP.nStims*txP.gabor.size_px);

        noiseCenters = [0:txP.gabor.size_px:(tkP.nStims-1)*txP.gabor.size_px; zeros(1, tkP.nStims)] + txP.gabor.size_px/2;
        baseRect = [0 0 txP.gabor.size_px txP.gabor.size_px];

        srcRects = CenterRectOnPointd(repmat(baseRect, [tkP.nStims,1])', noiseCenters(1,:), noiseCenters(2,:));

        minDist_px   = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.minDist_dva);
        minFixDist1 = txP.gabor.size_px*prm.fixDistFactor1;

        for b=1:tkP.nBlocks
            for i=1:nTrialsBuffered
        % (d) Matrizes com os centros dos estímulos (i.e., centros dos dstRects)
        %     e das cruzes de fixação
                [currFixCenter, currStimCenter] = getStimLocations(dpP.ellipseProps, tkP.nStims, minDist_px, prm.gridShape, prm.randomize);
                fixCenters(:, i, b)     = currFixCenter;
                stimCenters(:, :, i, b) = currStimCenter;
        % (e) Matriz com orientações tem os alvos adicionados
                orientation(randperm(tkP.nStims, nTs(b, i)), i, b) = targetOri(b);
            end
        end

        % (f) Distância mínima para considerar fixação em alvo pós-modificação
        minFixDist3 = txP.gabor.size_px*prm.fixDistFactor3;

        auxFixQueue = zeros(1, tkP.nStims);

    %% 4) Início dos blocos e trials
        try
            Screen('TextFont', dpP.window, prm.textFont);
            keepGoingBlocks = true;
            savedNdone = [0 0];
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
        
                    while true
                        [~,~,keyCode] = KbCheck;
                        if keyCode(spaceKey)
                            KbReleaseWait;
                            repeatMessage = false;
                            break;
                        elseif keyCode(escapeKey)
                            KbReleaseWait;
                            [keepGoingBlocks, ~, ~, ~] = pauseHandle(keepGoingBlocks, [], [], [], wasRecording, tkP, dpP, drP, prm, 'block', debug, mode);
                            break
                        end
                        WaitSecs(.01);
                    end
                end

                blockOnset = Screen('Flip', dpP.window); %#ok<NASGU>
                if debug == 0 && mode >= 2 && keepGoingBlocks
                    Eyelink('Message',sprintf(prm.msg.on.blk{1}, b, tkP.nBlocks));
                end
                i = 1;
                % Reordeno os trials para que, caso reinicie o bloco, a
                % ordem seja diferente
                trialQueue = randperm(nTrialsBuffered);
                retryCount = zeros(1, nTrialsBuffered);

                keepGoingTrials = keepGoingBlocks;
                while i <= tkP.nTrials && keepGoingTrials
                    trialIdxUp = false;
                    seenStimsQueue{b, i} = [];
                    idx = trialQueue(i);
                    trialOrder(1, i, b) = idx;
                    restartTrial = false;

                    fprintf('ATENÇÃO: %d visitas até modificar\n', modTimes(b, idx))

        % (b) Obtém a mediana do vetor de tempos de fixação
                    medFixTime = median(tkP.fixQueue);
                    if mode == 1
                        maxTrialDur = prm.maxTrialDurFactorCursor*tkP.nStims*medFixTime;
                    else
                        maxTrialDur = prm.maxTrialDurFactor*tkP.nStims*medFixTime;
                    end
        % (c) Cria os retângulos de destino com base nas coordenadas dos
        %     centros dos estímulos
                    dstRects = CenterRectOnPointd(repmat(baseRect, [tkP.nStims,1])', stimCenters(1,:, idx, b), stimCenters(2,:, idx, b));
                    
                
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
                    auxNoiseMatrix = butterFilter(pinkNoise(txP.gabor.size_px, tkP.nStims*txP.gabor.size_px), txP.noiseLoCutFreq, txP.noiseHiCutFreq);
                    noiseMatrix(:,:) = prm.noiseAlpha*(prm.noiseContrast*rescale(auxNoiseMatrix, -prm.noiseAmplitude, prm.noiseAmplitude))+drP.grey;
                    noiseTex = Screen('MakeTexture', dpP.window, noiseMatrix, [], [], [], [], []);
        
        % (e) Cria as texturas dos ruídos com orientação, sem desenhá-las
                    for j=1:tkP.nStims
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
                                        [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, dpP, drP, prm, 'trial', debug, mode);
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
                                        [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, dpP, drP, prm, 'trial', debug, mode);
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
                                        [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, dpP, drP, prm, 'trial', debug, mode);
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
        %      antes de ocorrerem as modificações pré-sacádicas
                    counter = 0;                 % Número de estímulos visitados
                    flag = zeros(1, tkP.nStims); % Quantas vezes cada estímulo foi visitado
                    currStim = 0;
                    fixStartTime = NaN;
                    if keepGoingTrials && ~restartTrial
%                         checkFixOnset = trialOnset;
                        if debug > 0 && mode > 1
                            seenIdx = sort(randsample(tkP.nStims,modTimes(b, idx)))';
                            flag(seenIdx) = 1;
                            currIdx = randsample(seenIdx, 1);
                            KbPressWait;
                        else
                            % Esse runTrial está sujeito a duas condições, vide
                            % fim do while abaixo, e diz respeito ao número
                            % de estímulos vistos
                            runTrial = true;
                            while runTrial
                                WaitSecs(0.001);
                                tNow = GetSecs;
        
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
                                        if currStim ~= 0
                                            fixDur = tNow - fixStartTime;
                                
                                            if fixDur >= prm.minFixTime2
                                                if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.off.stm{1}); end
                                                % Apenas salva na fila a primeira
                                                % fixação em um estímulo
                                                if flag(currStim) == 0
                                                    counter = counter + 1;
                                                    fprintf('Iniciou a visita ao %d-ésimo estímulo\n', counter)
                                                    auxFixQueue(counter) = fixDur;
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
                                        fixStartTime = NaN;
                                    
                                    % Entra no else quando há estímulo fixado
                                    % (i.e., começa ou continua fixando)
                                    else
                                        % Se logo antes não havia estímulo fixado,
                                        % é o início da fixação
                                        if currStim == 0
                                            currStim = currIdx;

                                            if debug == 0 && mode >= 2
                                                if flag(currStim) == 0
                                                    Eyelink('Message',prm.msg.on.stm{1});
                                                else
                                                    Eyelink('Message',prm.msg.on.stm{2});
                                                end
                                            end

                                            %% IMPORTANTE: Se for iniciada uma fixação
                                            % num modTimes(b,i)-ésimo estímulo diferente, 
                                            % inicia incrementa o counter, pois isso encerrará
                                            % o while runTrial
                                            if flag(currStim) == 0 && counter == modTimes(b, idx) - 1
                                                counter = counter + 1;
                                                fprintf('Visita ao último %d-ésimo estímulo\n', counter)
                                            end
                                            fixStartTime = tNow;

                                        % Se antes havia um estímulo diferente,
                                        % deve terminar a fixação e começar outra
                                          % (improvável que seja usado, só com
                                          %  estímulos próximos E amostragem baixa)
                                        elseif currStim ~= currIdx
                                            fixDur = tNow - fixStartTime;
        
                                            if fixDur >= prm.minFixTime2
                                                % De novo, apenas salva na fila a 
                                                % primeira fixação em um estímulo
                                                if flag(currStim) == 0
                                                    counter = counter + 1;
                                                    auxFixQueue(counter) = fixDur;
                                                end
                                                flag(currStim) = flag(currStim) + 1;
                                            end
        
                                            currStim = currIdx;
                                            fixStartTime = tNow;
                                        end
                                        % Se ainda estiver fixando o mesmo estímulo
                                        % (i.e., currStim == stmIdx), não faz nada
                                    end
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
                                                [keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial] = pauseHandle(keepGoingBlocks, restartBlock, keepGoingTrials, restartTrial, wasRecording, tkP, dpP, drP, prm, 'trial', debug, mode);
                                                break;
                                            end
                                        end
                                    end
                                end
                                runTrial = (counter < modTimes(b, idx)) && ~restartTrial;
                            end
                            seenIdx = find(flag ~= 0);
                        end
                        notSeenIdx = find(flag == 0);
                    end
        
    
    %% 7) Início Fase 3: tela com ruído rosa
        % Apenas se a tela 
                    if mode == 1
                        Screen('Close', bg);
                        clear auxWin bg;
                    end
        
        % (l) Desenha ruído orientado em todos os estímulos menos o 
        %     atual -- as linhas comentadas servem para mudar apenas os 
        %     não vistos. Atrasa a apresentação para o estímulo durar
        %     medFixTime segundo antes de o ruído rosa substituí-lo
                    if ~restartTrial && keepGoingTrials
                        blinkIdx = setdiff(1:tkP.nStims, currIdx);

            % (i) Desenha os gratings, somando ambos os sinais
                        Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
            %             Screen('DrawTextures', dpP.window, oriPinkTex, srcRects(:,notSeenIdx), dstRects(:,notSeenIdx), orientation(notSeenIdx, idx, b), [], [], [], []);
                        Screen('DrawTextures', dpP.window, oriPinkTex, srcRects(:,blinkIdx), dstRects(:,blinkIdx), orientation(blinkIdx, idx, b), [], [], [], []);

                        if ~isempty(currIdx)
                            Screen('DrawTextures', dpP.window, noiseTex, srcRects(:, currIdx), dstRects(:, currIdx), orientation(currIdx, idx, b), [], [], [], []);
                            Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
                            Screen('DrawTexture', dpP.window, txP.gabor.tex, [], dstRects(:, currIdx), orientation(currIdx, idx, b), [], [], [], [], [], txP.gabor.props);
                        end
                    
            % (j) Desenha a abertura gaussiana
                        Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            %             Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects(:,notSeenIdx), orientation(notSeenIdx, idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
                        Screen('DrawTextures', dpP.window, txP.blob.tex, [], dstRects, orientation(:, idx, b), [], [], [0 0 0 1]', [], [], txP.blob.props);
                        Screen('Close', oriPinkTex); Screen('Close', noiseTex);
                        Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        
                        %% A tela não modificada deve ser exibida ao todo por 
                        % medFixTime - prm.pinkNoiseDur, mas devo descontar
                        % o tempo passado desde o início da fixação
                        auxT1 = (medFixTime - prm.pinkNoiseDur); auxT2 = (GetSecs - fixStartTime);
                        timeLeft = max(0, auxT1 - auxT2);
                        fprintf('Tempo permitido de fixação antes do rosa: %.4f\n', auxT1)
                        fprintf('Tempo transcorrido desde início da fixação: %.4f\n', auxT2)
                        fprintf('Tempo restante ate exibir ruído rosa: %.4f\n', timeLeft)
                
            % viii. Registra os tempos de início e fim da Fase 3
                        updateStimOnset = Screen('Flip', dpP.window, timeLeft);
                            if debug == 0 && mode >= 2, Eyelink('Message',prm.msg.on.P3); end
                        if mode == 1, img = Screen('GetImage', dpP.window); bg = Screen('MakeTexture', dpP.window, img); end
                        
                        % A fase 3 é encerrada se o estímulo fica tempo
                        % demais na tela ou quando o olho sai do último
                        % estímulo
                        while true 
                            tNow = GetSecs;
                            if tNow - updateStimOnset > prm.pinkNoiseDur
                                auxT1 = tNow - updateStimOnset;
                                fprintf('Fim ruído rosa por duração: %.4f\n', auxT1)
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
                                isCurrStim = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, idx, b)) <= minFixDist1;
                                if isCurrStim(currStim) == 0
                                    auxT1 = tNow - updateStimOnset;
                                    fprintf('Fim ruído rosa por dispersão: %.4f\n', auxT1)
                                    break;
                                end
                            end
                            WaitSecs(.0005);
                        end
                        updateStimOffset = Screen('Flip', dpP.window);
                        auxT1 = updateStimOffset - updateStimOnset;
                        fprintf('Tempo total de ruído rosa: %.4f\n', auxT1)

                        if debug == 0 && mode >= 2,  Eyelink('Message',prm.msg.off.P3); end
%                         updateStimOffset = Screen('Flip', dpP.window, updateStimOnset + prm.pinkNoiseDur);
            
            % ix. Verifica se há alguma fixação de duração mínima em estímulo 
            %     numa janela pós-modificação
                        fixOnset = updateStimOffset; currIdx = [];
            
                        if debug == 0 || mode == 1
                            maxDurReached = true;
                            tNow = GetSecs;
                            while tNow - updateStimOffset < prm.postModDur
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
                                    isCurrStim = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, idx, b)) <= minFixDist3;
                                    if isempty(currIdx)
                                        if any(isCurrStim)
                                            currIdx = find(isCurrStim, 1);
                                            fixOnset = tNow;
                                        end
                                    else
                                        if ~isCurrStim(currIdx)
                                            fixDur = tNow - fixOnset;
                                            if fixDur >= prm.minFixTime3
                                                seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currIdx; fixDur]];
                                                disp(['Trial ' num2str(idx) ': Visitou o alvo ' num2str(currIdx) ' pós-modificação']);
                                                maxDurReached = false;
                                                break
                                            end
                                        end
                                        
                                    end
                                end
                                WaitSecs(.001);
                                tNow = GetSecs;
                            end
                            if maxDurReached && ~isempty(currIdx)
                               seenStimsQueue{b, i} = [seenStimsQueue{b, i} [currIdx; prm.postModDur]];
                            end
                        end
            
                        if ~isempty(currIdx) && ismember(currIdx, notSeenIdx)
                            notSeenIdx(notSeenIdx == currIdx) = [];
                        end
    %% 8) Início Fase 4: reportar em quais posições havia alvos
            % x. Desenha placeholders para os estímulos aleatoriamente, obedecendo
            %    a identificação deles como visitados, atual e não visitados
                        if mode == 1
                            Screen('Close', bg); clear auxWin bg;
                        end
                        allTargets = nan(1,tkP.nStims); allColors2 = drP.allColors;
                        Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            
                        seenAux = datasample(setdiff(seenIdx, currIdx), min(length(seenIdx), nStimsToReport(1, idx, b)), 'Replace', false);
                        currAux = []; if nStimsToReport(2, idx, b) == 1, currAux = currIdx; end
                        notSeenAux = datasample(notSeenIdx, min(length(notSeenIdx), nStimsToReport(3, idx, b)), 'Replace', false);
                        % Tanto seenAux como notSeenAux devem ser linhas
                        % para concatenar 
                        orderToReportStimsCell = {seenAux, currAux, notSeenAux};
                        orderToReportStims = [orderToReportStimsCell{orderToReportSets(1, idx, b)} orderToReportStimsCell{orderToReportSets(2, idx, b)} orderToReportStimsCell{orderToReportSets(3, idx, b)}];
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
                        
                        if debug == 0 && mode >= 2
                            Eyelink('Message',prm.msg.off.P4);
                        end

                        feedback = (orientation(:,idx,b) == targetOri(b))' + allTargets;
                        feedback(rem(feedback,2) == 0) = 2; feedback(rem(feedback,2) == 1) = 0;
                        feedback = feedback/2;

                        trialFeedback{b, i} = [orderToReportStims; orderRemapped; feedback(orderToReportStims)];
                        
                        if mode < 3
                            foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW, feedback, drP.red, drP.green);
                        else
                            foragingFlip(dpP.window, stimCenters(:, :, idx, b), dstRects, orderToReportStims, txP.gabor.size_px, drP.allColors, allTargets, targetOri(b), drP.allPW);
                         end
                        WaitSecs(1);
                        
                        tkP.fixQueue = [tkP.fixQueue(modTimes(b, idx):end), auxFixQueue(1:(modTimes(b, idx)-1))];
                        trialOrder(2, i, b) = 1;

                        i = i + 1;
                        trialIdxUp = true;
                    else
                        Eyelink('Message',prm.msg.err.trl{1});
                        retryCount(trialQueue(i)) = retryCount(trialQueue(i)) + 1;

                        if retryCount(trialQueue(i)) > prm.maxRetries
                            Eyelink('Message',prm.msg.err.trl{2});
                            warning('Trial %d excede o máximo de repetições. Prosseguindo', trialQueue(i));
                            i = i + 1;
                            trialIdxUp = true;
                        else
                            % Faz com que o trial não terminado sempre vá para
                            % o fim da fila
                            trialQueue(i:end) = [trialQueue(i+1:end) trialQueue(i)];
                        end
                     end
    
            % xi. Interrompe o registro, pois ou a tela será atualizada ou
            %     acabaram os trials
                    if debug == 0 && mode >= 2 && keepGoingTrials
                        trialOffset = GetSecs; %#ok<NASGU>
                        Eyelink('Message',sprintf(prm.msg.off.trl{1}, i - trialIdxUp, tkP.nTrials));
                        Screen('Flip', dpP.window);
                        WaitSecs(0.1);
                        Eyelink('SetOfflineMode');
                        Eyelink('StopRecording');
                    else 
                        Screen('Flip', dpP.window);
                        WaitSecs(0.5);
                    end
                end
                if keepGoingTrials
                    if restartBlock
                        Eyelink('Message',prm.msg.err.blk);
                        restartBlock = false;
                    else
                        blockOffset = GetSecs; %#ok<NASGU>
                        Eyelink('Message',sprintf(prm.msg.off.blk{1}, b, tkP.nBlocks));
                        b = b+1;
                    end
                end

                if b == tkP.nBlocks+1 && keepGoingBlocks
                    Eyelink('Message',sprintf(prm.msg.off.ses{1}, suffix));
                    blocksCompleted = true;
                elseif b ~= tkP.nBlocks+1 && ~keepGoingBlocks
                    Eyelink('Message',sprintf(prm.msg.off.ses{2}, suffix));
                end
            end

            % Salva todos os arquivos relevantes
            savedNdone(2) = blocksCompleted;
            if debug == 0 && mode >= 3

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

                if Eyelink('IsConnected')
                    outPath = fullfile(prm.filePath, prm.outFolder);
                    if Eyelink('CheckRecording')
                        Eyelink('SetOfflineMode');
                        Eyelink('StopRecording');
                    end
                    Eyelink('CloseFile');
    
                    % Por algum motivo puxa várias vezes do EyeLink
                    ntimes = 1;
                    while ntimes <= 10
                        status = Eyelink('ReceiveFile');
                        if status > 0
                            break
                        end
                        ntimes = ntimes+1;
                    end
                    if status <= 0
                        warning('Arquivo do Eyelink não recebido!')
                    end
                    edfFile = [tkP.edfFile prm.edfExtension]; 
                    if exist(edfFile,'file')
                        movefile(edfFile, outPath);
                        disp('Arquivo movido para o diretório fornecido.')
                        savedNdone(1) = true;
                    else
                        warning('Arquivo do EyeLink não encontrado!')
                    end
                end
                outFile = fullfile(outPath, ['s' tkP.edfFile]);
                save(outFile,'tkP', 'dpP', 'drP', 'txP', 'prm', 'results');
            else
                disp('Se fosse sessão experimental, estaria salvo')
            end
        catch
            psychrethrow(psychlasterror);
            cleanup(dpP.window);
        end

    if mode < 3
        tkP.nBlocks = nBlocks;
        tkP.nTrials = nTrials;
        if debug == 0, HideCursor(dpP.window); end
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
        modTimes = randsample(1:nStims, nBlocks*nTrials, true, targetModTimePMF);
        modTimes = reshape(modTimes, nBlocks, nTrials);
%         modTimes = modTimes - 1;

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


function [fixCenter, stimCenters] = getStimLocations(ROIparams, nStims, minDist, gridShape, randomize)
% Se quiser randomizar, independentemente da relação entre nStims e 
% size(gridSize(:)), divide ao acaso os estímulos por célula do grid;
%
% Se nStims < size(gridSize(:)), inevitavelmente randomiza;
%
% Se nStims >= size(gridSize(:)) e não quiser randomizar, distribui
% uniformemente ao longo das células e o que resta randomiza
    if nargin < 4, gridShape = [2 4]; end
    if nargin < 5, randomize = true;  end
    
    gridSize = prod(gridShape);
    
    if randomize || nStims < gridSize
    %     meanStim = nStims/gridSize; auxRand = sqrt(meanStim)*randn(gridSize, 1) + meanStim;
        auxRand = rand(gridSize, 1);
        nStimsEachCell = floor(auxRand.*nStims./sum(auxRand));
        
        rest = nStims - sum(nStimsEachCell);
    else
        allocStimsEachCell = floor(nStims/gridSize);
        rest = nStims - allocStimsEachCell*gridSize;
        nStimsEachCell = ones(gridShape)*allocStimsEachCell;
    end
    for i = 1:rest
        idx = randi(gridSize);
        nStimsEachCell(idx) = nStimsEachCell(idx) + 1;
    end
    rectLim = [ROIparams(1)-ROIparams(3), ROIparams(1)+ROIparams(3);
               ROIparams(2)-ROIparams(4), ROIparams(2)+ROIparams(4)];
    cellXlims = linspace(rectLim(1,1), rectLim(1,2), gridShape(2)+1);
    cellYlims = linspace(rectLim(2,1), rectLim(2,2), gridShape(1)+1);

    nStimsEachCell = reshape(nStimsEachCell, gridShape);
    
    [X, Y] = meshgrid(cellXlims, cellYlims);
    
    stimCenters = [];
    for j=1:gridShape(2)
        for i=1:gridShape(1)
            a = [X(i,j), X(i+1,j+1)]; b = [Y(i,j), Y(i+1,j+1)];
            aux = disc_rejection_sampling(a, b, nStimsEachCell(i, j), minDist, 'ellipse', ROIparams, stimCenters);
                 
             stimCenters = [stimCenters aux]; %#ok<AGROW>
        end
    end
      fixCenter = disc_rejection_sampling(rectLim(1,:), rectLim(2,:), 1, minDist, 'ellipse', ROIparams,stimCenters);
end


function currTime = foragingFlip(win, centers, dstCoord, idxForOvals, s, colors, aT, tOri, pW, fb, wrongCol, rightCol)
    if nargin >= 10
        rightCount = sum(fb == 1);
        wrongCount = sum(fb == 0);
        colors(:, fb == 1) = repmat(rightCol', [1, rightCount]);
        colors(:, fb == 0) = repmat(wrongCol', [1, wrongCount]);
    end
    Screen('FrameOval', win, colors(:, idxForOvals), dstCoord(:, idxForOvals), pW(:, idxForOvals));

    targets = (aT == 1);
%     neutralStim = (aT == -1);
    nonTargets = (aT == 0);
    % Somo 90 ficar em relação à vertical, como está o resto do código...
    drawInclinedLines(win, centers(:,targets), s*.8, tOri + 90, colors(:,targets), pW(:, targets));
    drawDots(win, centers(:,nonTargets), s/5, colors(:,nonTargets));

    currTime = Screen('Flip', win);
    WaitSecs(0.001); 
end


function drawDots(windowPtr, centers, len, color)
    if ~isempty(centers)
        Screen('DrawDots', windowPtr, centers, len, color, [], 2);
    end
end


function drawInclinedLines(windowPtr, centers, len, ang, color, width)

    if ~isempty(centers)
        % Simplified version for uniform lines
        numLines = size(centers, 2);
        angRad = deg2rad(ang);
        hl = len/2;
        
        % Calculate offsets
        dx = hl .* cos(angRad);
        dy = hl .* sin(angRad);
        
        % Create interleaved start and end points
        xy = [centers(1,:) - dx; centers(2,:) - dy; 
              centers(1,:) + dx; centers(2,:) + dy];
        xy = reshape(xy, 2, 2 * numLines);
        
        Screen('DrawLines', windowPtr, xy, width, color(:,repelem(1:size(color,2), 2)));
    end
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


function [kgBlocks, rBlock, kgTrials, rTrial] = pauseHandle(kgBlocks, rBlock, kgTrials, rTrial, wasRecording, tkP, dpP, drP, prm, pauseMode, debug, foragingMode)
    pauseStart = GetSecs;
    Eyelink('Message',prm.msg.on.pse);
    Eyelink('StopRecording');
    Eyelink('SetOfflineMode');
    decision = pauseScreen(tkP, dpP, drP, prm, pauseMode);
    
    if strcmp(pauseMode, 'trial')
        if strcmp(decision, 'gotoMenu')
            Eyelink('Message',prm.msg.pse{1});
            kgTrials = false;
            kgBlocks = false;
        elseif strcmp(decision, 'restartBlock')
            Eyelink('Message',prm.msg.pse{2});
            rBlock = true;
            kgTrials = false;
        elseif strcmp(decision, 'recalibrate')
          Eyelink('Message', prm.msg.pse{3});
            if debug == 0 && foragingMode >= 2
                EyelinkDoTrackerSetup(tkP.el);
            else
                disp('Recalibragem solicitada')
            end
            rTrial = true;
        elseif strcmp(decision, 'resume')
            Eyelink('Message',prm.msg.pse{4});
            rTrial = true;
        end
        pauseEnd = GetSecs;
        pauseDur = pauseEnd - pauseStart;
        if kgTrials
            if pauseDur > prm.pauseMaxDur
                fakeLoadingScreen(tkP, dpP, drP, prm, 'trial', txP, targetOri(b));
                Screen('Flip',dpP.window);
            end
            WaitSecs(.01);
            if wasRecording
                Eyelink('StartRecording'); 
            end
        end
    elseif strcmp(pauseMode, 'block')
        if strcmp(decision, 'gotoMenu')
            Eyelink('Message', prm.msg.pse{1});
            kgBlocks = false;
        elseif strcmp(decision, 'recalibrate')
            Eyelink('Message',prm.msg.pse{3});
            if debug == 0 && mode >= 2
                EyelinkDoTrackerSetup(tkP.el);
            else
                disp('Recalibragem solicitada')
            end
        elseif strcmp(decision, 'resume')
            Eyelink('Message',prm.msg.pse{4});
        end
        if kgBlocks
            if wasRecording
                Eyelink('StartRecording'); 
            end
        end
    end
    Eyelink('Message',prm.msg.off.pse);
end


function decision = pauseScreen(tkP, dpP, drP, prm, mode)
    decision = '';

    leaveOptions = {'quit', 'gotoMenu'};
    parentDir = prm.filePath;
    leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; spaceKey = tkP.keys{3};
    cX = dpP.winCenter(1); cY = dpP.winCenter(2);


    titleText = 'PAUSADO'; titleMargin = cX*prm.titleMarginFactor;
    %% Tipos de botões
    switch mode
        case 'menu'
            options     = {'quit', 'resume'};
            optionsName = {'Sair', 'Prosseguir'};

        case 'trial'
            options = {'gotoMenu', 'restartBlock', 'recalibrate', 'resume'};
            optionsName = {'Menu', 'Reiniciar Bloco', 'Recalibrar', 'Prosseguir'};

        case 'block'
            options = {'gotoMenu', 'recalibrate', 'resume'};
            optionsName = {'Menu', 'Recalibrar', 'Prosseguir'};
    end

    icons    = cellfun(@(x) fullfile(parentDir, prm.imgFolder, [x prm.imgExtension]), options, 'UniformOutput', false);
    iconsImg = cellfun(@(x) alphaRead(x), icons, 'UniformOutput', false);
    imgSizes = cellfun(@(x) size(x), iconsImg, 'UniformOutput', false);
    auxImg   = cellfun(@(x) drP.black*ones(x(1), x(2), 4), imgSizes, 'UniformOutput', false);
    iconsImg = cellfun(@(x, a) addAlpha(x,a), auxImg, iconsImg, 'UniformOutput', false);
    iconsTex = cellfun(@(x) Screen('MakeTexture', dpP.window, x), iconsImg);
    clear icons imgSizes auxImg iconsImg;

    L = numel(options);
    if (L+1)/2 == floor((L+1)/2)
        firstSelectL = (L+1)/2;
        firstSelectR = (L+1)/2;
    else
        firstSelectL = floor((L+1)/2);
        firstSelectR = ceil((L+1)/2);
    end
    selected = 0;

    %% Retângulos dos botões
    gap  = 80;
    btnW = min(300, (dpP.winRect(3)-(L+1)*gap)/L); btnH = 300;
    totalW = L*btnW + (L-1)*gap;
    leftX0 = cX - totalW/2;
    leftX  = leftX0 + (0:L-1) * (btnW + gap);

    btnRects = zeros(4, L);
    btnRects(1,:) = leftX; btnRects(2,:) = cY - btnH/2;
    btnRects(3,:) = leftX + btnW; btnRects(4,:) = cY + btnH/2;
    
    prm.upFrac = 2/3; SizeFactor = .6;
    RectSplit = btnRects(2,:)+btnH*prm.upFrac;
    upRects = btnRects; upRects(4,:) = RectSplit; downRects = btnRects; downRects(2,:) = RectSplit;
    uRH = RectHeight(upRects(:,1)'); uRW = RectWidth(upRects(:,1)');
    uRHW = min(uRH, uRW);
    for i=1:L
        [cXaux, cYaux] = RectCenter(upRects(:,i));
        upRects(:,i) = CenterRectOnPoint([0 0, uRHW, uRHW]*SizeFactor, cXaux, cYaux);
    end
    clear cXaux cYaux

    keyWasDown = false;
    lastActionTime = 0; lastRepeatTime = -Inf;
    doAction = false;

    KbReleaseWait;

    while true
        Screen('TextStyle', dpP.window, 0);
        Screen('TextSize', dpP.window, prm.textSizeTitle);
        Screen('DrawText', dpP.window, titleText, titleMargin, titleMargin,  drP.blackGrey);

        Screen('TextSize', dpP.window, prm.textSizeNormalish);

        %% Desenha os retângulos e textos e texturas
        Screen('TextStyle', dpP.window, 1);
        Screen('FillRect', dpP.window, drP.blue, btnRects);
        Screen('DrawTextures', dpP.window, iconsTex, [], upRects);

        for i=1:L
            DrawFormattedText(dpP.window, optionsName{i}, 'center', 'center',  drP.black, [], [], [], [], [], downRects(:,i)');
        end

        % Se houver algum retângulo selecionado, adiciona contorno
        if selected ~= 0
            Screen('FrameRect', dpP.window,  drP.darkBlue, btnRects(:, selected), prm.pW2);
        end
        Screen('Flip', dpP.window);

        [keyIsDown, ~, keyCode] = KbCheck;

        if keyIsDown
            if ~keyWasDown
                lastActionTime = GetSecs;
                doAction = true;
            elseif GetSecs - lastActionTime > prm.repeatDelay
                if GetSecs - lastRepeatTime > prm.repeatRate
                    doAction = true;
                    lastRepeatTime = GetSecs;
                end
            end
            if doAction
                Screen('TextStyle', dpP.window, 0);
                doAction = false;
                if keyCode(leftKey)
                    if selected == 0
                        selected = firstSelectL;
                    else
                        selected = max(1, selected - 1);
                    end
    
                elseif keyCode(rightKey)
                    if selected == 0
                        selected = firstSelectR;
                    else
                        selected = min(L, selected + 1);
                    end
    
                elseif keyCode(spaceKey)
                    if selected ~= 0
                        decision = options{selected};
                    end
                    KbReleaseWait
                    if ismember(decision, leaveOptions)
                        KbReleaseWait;
                        confirmed = confirmationScreen(dpP, drP, prm, tkP.keys);
                        if confirmed
                            return;
                        else
                            decision = '';
                        end
                    end
                    if ~isempty(decision)
                        for i=1:L, Screen('Close', iconsTex(i)); end
                        return;
                    end
                end
            end
        else
            lastRepeatTime = -Inf;
        end
        keyWasDown = keyIsDown;
        WaitSecs(0.01);
    end
end


function alpha = alphaRead(x)
    [~,~,alpha] = imread(x);
end


function x = addAlpha(x, a)
    x(:,:,4) = a;
end


function confirmed = confirmationScreen(dpP, drP, prm, keys)

    options  = {'Sim', 'Cancelar'};
    selected = 2;
    blueRect = CenterRectOnPointd([0,0,1.2*dpP.winCenter(1),1.2*dpP.winCenter(2)], dpP.winCenter(1),dpP.winCenter(2));

    %% Posições e tamanhos
    msgY1     = dpP.winCenter(2) * 0.75;
    msgY2     = dpP.winCenter(2) * 0.85;
    btnY      = dpP.winCenter(2) * 1.35;
    btnW      = dpP.winRect(3) * 0.15;
    btnH      = dpP.winRect(4) * 0.08;
    btnGap    = dpP.winRect(3) * 0.1;

    cX = dpP.winCenter(1);

    % Retângulos simétricos em relação ao eixo vertical
    btnRects(:, :, 1) = CenterRectOnPoint([0 0 btnW btnH], cX - btnGap, btnY);
    btnRects(:, :, 2) = CenterRectOnPoint([0 0 btnW btnH], cX + btnGap, btnY);

    KbReleaseWait;
    while true
        
        Screen('FillRect', dpP.window, drP.blue, blueRect);
        Screen('TextSize', dpP.window, prm.textSizeBigger);
        DrawFormattedText(dpP.window, 'Tem certeza que deseja sair?\n', 'center', msgY1, drP.black);
        Screen('TextSize', dpP.window, prm.textSizeBig);
        DrawFormattedText(dpP.window, 'Essa decisão não poderá ser desfeita', 'center', msgY2, drP.blackGrey);

        for i = 1:2
            if i == selected
                frameColor = drP.red;
                frameWidth = 4;
            else
                frameColor = drP.black;
                frameWidth = 2;
            end

            Screen('FrameRect', dpP.window, frameColor, btnRects(:,:,i), frameWidth);

            DrawFormattedText(dpP.window, upper(options{i}), ...
                'center', 'center', drP.black, [], [], [], [], [], btnRects(:,:,i));
        end

        Screen('Flip', dpP.window);

        % Teclado
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(keys{1})     % Seta para esquerda
                selected = 1;

            elseif keyCode(keys{2}) % Seta para direita
                selected = 2;

            elseif keyCode(keys{3}) % Espaço (confirmar)
                confirmed = strcmp(options{selected}, 'Sim');
                Screen('TextSize', dpP.window, prm.textSizeNormal);
                KbReleaseWait;
                return;
            end
            KbReleaseWait;
        end
    end
end


function fakeLoadingScreen(tkP, dpP, drP, prm, mode, txP, ori)
    % Tela de 'Carregando' falsa com barra e dicas

    if nargin < 5, mode = 'opening'; end

    leftKey = tkP.keys{1}; rightKey = tkP.keys{2}; escapeKey = tkP.keys{4};

    parentDir = fileparts(mfilename('fullpath'));
    fileLogo1 = 'logo_bBg_wSb.png'; fileLogo2 = 'logo_gBg_wSb.png'; fileRoundRect = '1by2_rect.png';

    pathLogo1 = fullfile(parentDir, prm.imgFolder, fileLogo1);
    pathLogo2 = fullfile(parentDir, prm.imgFolder, fileLogo2);

    pathRectRound = fullfile(parentDir, prm.imgFolder, fileRoundRect);
    clear parentDir fileRoundRect

    if strcmp(mode, 'opening')
        pathLogo = pathLogo1; f = 0;
        bgColor   = drP.black;
        totalTime = 15;
        hints     = {
                     'Ganhe brinde extra ao indicar outras pessoas\n para participarem de nossos experimentos!'
                     'Faça parte da lista de sujeitos experimentais\n do laboratório!'
                     'Evite distrações enquanto realiza a tarefa.'
                     'Se o experimento tiver mais de uma sessão, \nnão esqueça de comparecer em todas elas!'
                    };

        hintChangeInterval = 6; % Tempo durante o qual uma dica é exibida
        showBar     = true;
        barYOffset  = 135;
        hintYOffset = 75;       % Quantos pixels abaixo da barra aparecem as dicas
    elseif ismember(mode, {'start', 'trial'})
        pathLogo = pathLogo2; f = 1;
        bgColor   = drP.grey;
        totalTime = 8;
        hints     = {
                     'Memorize apenas onde estão os alvos!'
                     'Faça pausas apenas quando o tipo de alvo mudar.'
                     'Evite se mexer durante o experimento.'
                     'Tente se manter concentrado e relaxado.'
                     'Priorize o desempenho, e não a velocidade.'
                     'Não esqueça de olhar para a cruz no início\nde cada tentativa!'
                    };
        hintChangeInterval = 4.5;
        showBar     = false;
        barYOffset  = 105;
        hintYOffset = 105;
    end
    
    %% Parâmetros
    nBlocks          = 20;      % Número de blocos na barra de progresso
    nSteps           = 20;      % São dados tantos passos quanto os blocos
    barHeight        = 20;      % Altura da barra de progresso em pixels
    blockGap         = 6;       % Espaço entre blocos em pixels
    imageYOffsetFrac = 0.15;    % Quão para cima ficará a imagem em relação ao centro
    creditXOffset    = 50;
    creditYOffset    = 25;
    arrowsXOffset    = 10;

    blockImg = imread(pathRectRound);
    
    white = repmat(drP.white, 1, 3);
    
    Screen('TextFont', dpP.window, prm.textFont);
    Screen('TextColor', dpP.window, drP.white);
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    cX = dpP.winCenter(1); cY = dpP.winCenter(2);
    
    %% Carrega a imagem do logo e ajusta seu tamanho (no máximo 2/3 da tela)
    
    logoImg = imread(pathLogo);
    % Cria canal alfa para remover o fundo (f é parâmetro ajustado conforme 
    % a cor do fundo)
    if size(logoImg,3) == 3
        meanLogo = mean(logoImg, 3);
        alpha = uint8( 255*(meanLogo-f*min(meanLogo(:)))/( max(meanLogo(:)) - min(meanLogo(:)) ) );
        logoImg(:,:,4) = alpha;
    end
    logoTex = Screen('MakeTexture', dpP.window, logoImg);
    clear pathLogo1 pathLogo2 pathLogo pathRectRound
    
    [imgH, imgW, ~] = size(logoImg);
    maxImgW = dpP.winRect(3) * 1/3; maxImgH = dpP.winRect(4) * 1/3;
    
    scaleFactor = min([1, maxImgW/imgW, maxImgH/imgH]);
    scaledW = imgW * scaleFactor; scaledH = imgH * scaleFactor;
    
    imgRect = CenterRectOnPointd([0 0 scaledW scaledH], cX, cY - dpP.winRect(4) * imageYOffsetFrac);
    clear logoImg imgH imgW maxImgW maxImgH scaleFactor scaledW scaledH

    
    if strcmp(mode, 'trial')
        totalTime = totalTime*(3/8);
        Screen('Close', logoTex);
        logoTex = -1;
    end

    %% Retângulo dos créditos
    creditX = imgRect(1) - creditXOffset; creditY = imgRect(2) - creditYOffset;
    if logoTex == -1
        creditText = 'ORIENTAÇÃO DO ALVO:';
    else
        creditText = 'UM EXPERIMENTO DE';
    end
    
    %% Carrega a imagem do retângulo arredondado como textura
        % Se não tiver canal alfa, cria de modo que o preto seja
        % transparente (para usar a cor do fundo)
    if size(blockImg,3) == 3
        alpha = uint8(mean(blockImg, 3));
        blockImg(:,:,4) = alpha;
    end
    
    blockTex = Screen('MakeTexture', dpP.window, blockImg);
    
    %% Barra de progresso
    % Tamanho original da imagem do bloco
    [bH, bW, ~] = size(blockImg);
    blockAspect = bW / bH;
    
    % A altura da imagem é dada com base na altura da barra
    blockHeight = barHeight;
    blockWidth  = blockHeight * blockAspect;
    
    % A largura da barra é dada com base no tamanho dos blocos e no espaço
    % entre eles
    totalBlocksWidth = nBlocks * blockWidth + (nBlocks-1) * blockGap;
    barX = cX - totalBlocksWidth / 2;
    barY = cY + barYOffset;
    
    % Contorno da barra
    barOutline = [barX - blockGap, barY - blockGap, barX + totalBlocksWidth + blockGap, barY + barHeight + blockGap];
    [barCX, barCY] = RectCenter(barOutline);
    
    % Retângulos em que são apresentadas as imagens dos blocos
    blockRects = zeros(4, nBlocks);
    for i = 1:nBlocks
        x1 = barX + (i-1) * (blockWidth + blockGap);
        blockRects(:,i) = [x1 barY x1+blockWidth barY+blockHeight];
    end
    clear blockImg bH bW blockAspect blockHeight blockWidth totalBlocksWidth
    clear barHeight barX blockGap imageYOffsetFrac creditXOffset creditYOffset 
    
    % As durações de cada passo são sorteadas uniformemente entre 0 e 1, de
    % modo que a soma coincida com a duração total
    stepDurations = rand(1, nSteps);
    stepDurations = stepDurations / sum(stepDurations) * totalTime;
    stepTimes = cumsum(stepDurations);
    
    % Aumenta um pouquinho o tempo total para ter algum instante com o 10/10
    totalTime = totalTime + prm.repeatDelay;
    blocksPerStep = nBlocks / nSteps;
    
    %% Dicas
    Screen('TextSize', dpP.window, prm.textSizeNormal);
    hintY = barY + hintYOffset;
    maxW = 0; maxH = 0;
    for i = 1:numel(hints)
        [~, ~, hintsRect] = DrawFormattedText(dpP.window, hints{i});
        currW = RectWidth(hintsRect); currH = RectHeight(hintsRect);
        if currW > maxW, maxW = currW; end
        if currH > maxH, maxH = currH; end
    end
    hintsRect = CenterRectOnPointd([0, 0, maxW, maxH], cX, hintY);
    clear stepDurations maxW maxH currW currH

    %% Carregando alternativo
    nDots = 4; amp = 8; freq = 2; k = 3;
    phaseStep = 2*pi/nDots; phase = (0:(nDots-1))*phaseStep;
    adjSin = @(x) amp .* sin(2*pi*freq*x + phase) .* (mod(floor(freq*x), k) == 0) .*(0.5 - 0.5*cos(2*pi*mod(freq*x,1)));
    radius = 3; distance = 3;
    loadingText = 'CARREGANDO';
    Screen('TextSize', dpP.window, prm.textSizeBig);
    boundsLoading = Screen('TextBounds', dpP.window, loadingText);
    loadingRect = CenterRectOnPointd(boundsLoading, barCX, barCY);
    loadingWidth  = RectWidth(boundsLoading);
    boundsAux = Screen('TextBounds', dpP.window, 'A');
    Screen('TextSize', dpP.window, prm.textSizeNormal);

    dotsX = cX + loadingWidth/2 + cumsum(repmat(distance+radius, 1, nDots));
    dotsY = repmat(loadingRect(2) + RectHeight(boundsAux), 1, nDots);

    %% Setas
    lArrowText = '<'; rArrowText = '>';
    arrowsRect = Screen('TextBounds', dpP.window, lArrowText);
    arrowsWidth = RectWidth(arrowsRect);
    larrowRect = CenterRectOnPointd(arrowsRect, hintsRect(1)-arrowsWidth-arrowsXOffset, hintY);
    rarrowRect = CenterRectOnPointd(arrowsRect, hintsRect(3)+arrowsWidth+arrowsXOffset, hintY);
    clear boundsLoading loadingWidth boundsAux nDots arrowsXOffset arrowsRect
    
    %% Loop principal
    hintIdx = 1;
    hints = hints(randperm(numel(hints)));
    currentHint = hints{mod(hintIdx, numel(hints))+1};
    hintPhase = 'steady';   % 'fadein', 'steady', 'fadeout'
    steadyStart = 0;
    hintAlpha = 1;
    
    if showBar, filledBlocks = 0; else, filledBlocks = nBlocks; end
    startTime = GetSecs;
    stepIndex = 0;
    
    doAction = false; selected = 0; keyWasDown = 0;
    Screen('TextSize', dpP.window, prm.textSizeSmall);
    lastActionTime = GetSecs;
    while true
        tNow = GetSecs;
        elapsedTime = tNow - startTime;
    
        % A tela é interrompida apenas se excedeu o limite de tempo e se
        % todos os blocos foram preenchidos
        if elapsedTime >= totalTime && filledBlocks == nBlocks
            break;
        end
    
        if showBar
            % Atualiza a quantidade de blocos preenchidos
            if stepIndex < nSteps && elapsedTime >= stepTimes(stepIndex + 1)
                stepIndex = stepIndex + 1;
            end
            filledBlocks = round(stepIndex * blocksPerStep);
        end

        %% Setas de seleção das dicas
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if ~keyWasDown
                lastActionTime = GetSecs;
                doAction = true;
            elseif GetSecs - lastActionTime > prm.repeatDelay
                if GetSecs - lastRepeatTime > prm.repeatRate
                    doAction = true;
                    lastRepeatTime = GetSecs;
                end
            end
            if doAction
                doAction = false;
                if keyCode(leftKey)
                    selected = -1;
                    hintPhase = 'fadeout'; fadeStart = elapsedTime;
                elseif keyCode(rightKey)
                    selected = 1;
                    hintPhase = 'fadeout'; fadeStart = elapsedTime;
                elseif keyCode(escapeKey)
                    break;
                end
            end
        end
    
        %% Atualiza a dica exibida
        % Se a dica ficou por tempo suficiente, ela deve sair e vir outra
        switch hintPhase
            case 'fadein'
                t = (elapsedTime - fadeStart) / prm.fadeInDur1;
                hintAlpha = min(1, t);
                % Se acabou a fase de transição (normalizada), muda
                % para estado steady
                if t >= 1
                    hintPhase = 'steady';
                    steadyStart = elapsedTime;
                end
        
            case 'steady'
                hintAlpha = 1;
        
                if elapsedTime - steadyStart >= hintChangeInterval - 2*prm.fadeInDur1
                    hintPhase = 'fadeout';
                    fadeStart = elapsedTime;
                end
        
            case 'fadeout'
                t = (elapsedTime - fadeStart) / prm.fadeInDur1;
                hintAlpha = max(0, (1 - t));
        
                if t >= 1
                    hintPhase = 'fadein';
                    fadeStart = elapsedTime;
                    if selected ~= 0
                        hintIdx = hintIdx + selected;
                        selected = 0;
                    else
                        hintIdx = hintIdx + 1;
                    end
                    currentHint = hints{mod(hintIdx, numel(hints))+1};
                end
        end

        Screen('FillRect', dpP.window, bgColor);
    
        %% Desenha o logo e o texto logo acima
        
        if logoTex == -1
            Screen('TextSize', dpP.window, prm.textSizeNormal);
            Screen('BlendFunction', dpP.window, GL_ONE, GL_ZERO);
            Screen('DrawTexture', dpP.window, txP.exampleNoise.tex, [], imgRect, ori, [], [], [], []);
            Screen('BlendFunction', dpP.window, GL_ONE, GL_ONE);
            Screen('DrawTexture', dpP.window, txP.exampleGabor.tex, [], imgRect, ori, [], [], [], [], [], txP.exampleGabor.props);
            Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTexture', dpP.window, txP.exampleBlob.tex, [], imgRect, ori, [], [], [0 0 0 1]', [], [], txP.exampleBlob.props); 
            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        else
            Screen('DrawTexture', dpP.window, logoTex, [], imgRect);
        end
        DrawFormattedText(dpP.window, creditText, creditX, creditY, drP.white);
    
        %% Desenha a barra de progresso
        Screen('TextSize', dpP.window, prm.textSizeNormal);
        if showBar
            Screen('FrameRect', dpP.window, drP.white, barOutline, 2);
            if filledBlocks > 0
                Screen('DrawTextures', dpP.window, blockTex, [], blockRects(:,1:filledBlocks));
            end
        
            % Texto das etapas concluídas
            progressText = sprintf('%d/10%', floor(10 * filledBlocks / nBlocks));
            DrawFormattedText(dpP.window, progressText, 'center', barY - 25, drP.white);
        else
            Screen('TextSize', dpP.window, prm.textSizeBig);
            DrawFormattedText(dpP.window, loadingText, 'center', 'center', drP.white, [], [], [], [], [], loadingRect);
            Screen('DrawDots', dpP.window, [dotsX; dotsY+adjSin(GetSecs - startTime)], radius, drP.white, []);
            Screen('TextSize', dpP.window, prm.textSizeNormal);
        end

        %% Desenha/escreve a dica
        DrawFormattedText(dpP.window, currentHint, 'center', 'center', [white hintAlpha], [], [], [], [], [], hintsRect);
%         Screen('FrameRect', dpP.window, [], hintsRect);

        %% Desenha as setas
        Screen('TextSize', dpP.window, prm.textSizeSmall);
        Screen('TextStyle', dpP.window, max(-1*selected, 0));
        DrawFormattedText(dpP.window, lArrowText, 'center', 'center', drP.white, [], [], [], [], [], larrowRect);
        Screen('TextStyle', dpP.window, max(selected, 0));
        DrawFormattedText(dpP.window, rArrowText, 'center', 'center', drP.white, [], [], [], [], [], rarrowRect);
        Screen('TextStyle', dpP.window, 0);
%         Screen('TextSize', dpP.window, prm.textSizeNormal);
    
        Screen('Flip', dpP.window);
    end
    
    %% Faz uma transição lenta à cor original
    if logoTex ~= -1, Screen('Close', logoTex); end
    Screen('Close', blockTex);

    img = Screen('GetImage', dpP.window);
    frozenTex = Screen('MakeTexture', dpP.window, img);

    fadeStart = GetSecs;
    while true
        elapsed = GetSecs - fadeStart;
        t = min(elapsed / prm.fadeInDur1, 1);
        texAlpha = 1-t;
    
        Screen('FillRect', dpP.window, drP.grey);
        Screen('DrawTexture', dpP.window, frozenTex, [], [], [], [], texAlpha);
    
        Screen('Flip', dpP.window);
    
        if t >= 1
            break;
        end
    end
    
    Screen('Close', frozenTex);
    WaitSecs(prm.fadeInDelay1);
end

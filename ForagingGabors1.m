function ForagingGabors1(nTrials, nStims, nBlocks, nTargets, options)
        arguments
            nTrials {mustBeNumeric}   = 20;
            nStims {mustBeNumeric}    = 8;
            nBlocks {mustBeNumeric}   = 15;
            nTargets {mustBeNumeric}  = ceil(nStims/2);
            options.mode string       = 'experiment'
        end
        debug = 0;
        if strcmp(options.mode, 'debug'), debug = 1; elseif strcmp(options.mode, 'debugTV'), debug = 2; end
        disp(debug)
        if debug > 0
            nTrials  = 4;
            nStims   = 8;
            nBlocks  = 2;
            nTargets = ceil(nStims/2);
        end
        
        % Os alvos são apresentados com uma distribuição triangular
        % discreta ao redor da média nTargets, de modo que a quantidade de
        % estímulos não vistos terá em média uma metade composta por alvos
        if nTargets > 1 && nTargets < nStims
            numTargetsPMF = [.25 .50 .25];
            targetsSpace  = nTargets + [-1 0 1];
            As = randsample(targetsSpace, nBlocks*nTrials, true, numTargetsPMF);
            As = reshape(As, [nBlocks nTrials]);
        else
            As = nTargets*ones(nBlocks, nTrials);
        end
        cleanup
    
    
%% 1) Parâmetros das condições experimentais
    % (a) Distância da tela, em cm
        if debug == 2, screenDist = 200; else, screenDist = 57; end

    % (b) Tamanhos e frequências, em dva e cpd
        gaborSize_dva    = 1.5;
        gaborFreq_cpd    = 11; % usa 11
        noiseCutFreq_cpd = 5; % em cpd

        noiseLoCutFreq_cpd = 1/gaborSize_dva;  % em cpd
        noiseHiCutFreq_cpd = max(5,2*noiseLoCutFreq_cpd);  % em cpd

    % (c) Disposição e orientação básicas dos estímulos
        minDist_dva = max(3, gaborSize_dva);% OBS: ajustar com ellipseToScreenRatio 
        ellipseToScreenRatio = [3/4 3/4];   % Proporção da elipse em relação ao
                                              % tamanho da tela em cada dimensão
        gridShape   = [2 4];                % #linhas x #colunas
        randomize = false;
        allOri = [0 45 90 135];
        allOriMap  = containers.Map(allOri,1:length(allOri));
        allOriName = {'vertical', 'diagonal crescente', 'horizontal', 'diagonal decrescente'};
        % Uma orientação-alvo para cada bloco
        targetOri = randsample(allOri, nBlocks, true);

    % (d) Parâmetros de fixação
        % Tempo mínimo de fixação na cruz inicial
        minFixTime1 = .5;  % em s
        fixDistFactor1 = 1.2;
        %% Tempo médio (?) durante o qual o participante fixa nos estímulos
        % durante a busca. Será usado para garantir que em alguns trials se
        % observe o efeito pré-sacádico. O ideal é que seja distribuição 
        % com média e variância adapatadas para cada participante
        minFixTime2 = .14; % em s
        %% Tempo de fixação mínimo no próximo estímulo durante a tela cinza 
        minFixTime3 = minFixTime2/2; % em s
        % Maior tolerância a erro de fixação, já que o estímulo é removido
        postModDur  = .3; % em s
        fixDistFactor3 = 1.5;
        
    % (e) Distribuição do instante de modificação dos demais estímulos
        modeDistr = floor(3*nStims/4); rng('shuffle');
        %% Arrumar a função robust_beta_pmf
        targetModTimePMF = robust_beta_pmf(nStims, modeDistr, 'peakness', 1.5);
        auxNs = randsample(1:nStims, nBlocks*nTrials, true, targetModTimePMF);
        auxNs = reshape(auxNs, nBlocks, nTrials);
        Ns = auxNs - 1;

    % (f) Parâmetros do ruído rosa com orientação
        %% Arrumar parâmetro
        aSigma = 20; rSigma2 = .76;
        % Duração do estímulo de ruído rosa, em segundos
        pinkNoiseDur = .15;

    % (g) Quantidade de estímulos cuja orientação deve ser reportada
        propTrialsPSA = .5;     % Em metade dos trials espera-se que a condição
                                % permita avaliar atenção pré-sacádica
        minToReport = 2; maxToReport = 3;
        nStimsToReport = round((maxToReport - minToReport)*rand(nBlocks,nTrials) + minToReport);
        reportPostModFix = rand(nBlocks, nTrials) <= propTrialsPSA;
        nStimsToReport = nStimsToReport - reportPostModFix; % Guarda a quantidade de estímulos
                                                            % a serem reportados que não são 
                                                            % alvos da fixação 
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
%         correctnSTR1 = auxNs >= nStimsToReport(1,:);
%         correctnSTR3 = nStims - auxNs <= nStimsToReport(3,:);
            % Em cada trial, o sujeito deve reportar dos estímulos vistos
            % anteriormente, atual ou não vistos em ordem aleatória
        [~, orderToReportSets] = sort(rand(3, nTrials, nBlocks), 1);
    
    
%% 2) Inicializa PTB e EyelinkToolBox
    % i. Caixa de diálogo
        prompt = {'Sujeito', 'Número da sessão', 'Olho dominante (E ou D)'};
        dlg_title = 'Informações da sessão experimental';
        def = {'00', '01', 'R'};
        options.Resize = 'on';
        answer = inputdlg(prompt, dlg_title, 1, def, options);
        domEye = answer{3};
        if  isempty(answer), fprintf('Session cancelled by user\n'); cleanup; return; end
        edfFile = [answer{1} '_' answer{2}];

        % Interrompe se nome muito longo
        if length(edfFile) > 8, fprintf('Filename needs to be no more than 8 characters long (letters, numbers and underscores only)\n'); cleanup; return; end
    
    
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
    

    % (b) Escolhe a tela em que haverá o desenho e define algumas cores
        screenNumber = max(Screen('Screens'));
        white = WhiteIndex(screenNumber);
        grey = white / 2;
        [window, ~] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
                                [], [],  kPsychNeed32BPCFloat);

    % (c) Obtém propriedades da tela
        [monitorW_mm, ~] = Screen('DisplaySize', screenNumber); % Tamanho     da tela em mm
        screenRes = Screen('Resolution', screenNumber);         % Resolução   da tela em px
        ifi = Screen('GetFlipInterval', window);

    % (d) Obtém o centro e os semieixos da elipse
        ellipseProps = [screenRes.width/2, screenRes.height/2];
        ellipseProps = [ellipseProps ellipseProps(1)*ellipseToScreenRatio(1) ellipseProps(2)*ellipseToScreenRatio(2)];
        Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % (e) Define outras cores de interesse e espessura das canetas para
    %     desenhar algumas das figuras
        midGrey = (white+grey)/2;
        orange = [250, 140, 15]'/255;
        pW1 = 3; pW2 = 5;

        allColors = white*ones(3, nStims);
        allPW     = pW1*ones(1,nStims);

    % i. Inicializa Eyelink
        dummymode = 0;
        if debug > 0, dummymode = 1; disp('Debug: O EyeLink não será inicializado'); end
        EyelinkInit(dummymode); % Initialize EyeLink connection
%         status = Eyelink('IsConnected');
%         if status < 1, dummymode = 1; end

    % ii. Abre o arquivo .edf
        failOpen = Eyelink('OpenFile', edfFile);
        % Interrompe se arquivo não abrir
        if failOpen ~= 0, fprintf('Cannot create EDF file %s', edfFile); cleanup; return; end

    % iii. Obtém versão do EyeLink
%         ELsoftwareVersion = 0; [ver, versionstring] = Eyelink('GetTrackerVersion');
%         if ver ~=0, [~, vnumcell] = regexp(versionstring,'.*?(\d)\.\d*?','Match','Tokens'); ELsoftwareVersion = str2double(vnumcell{1}{1}); end

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
    
        % Put EyeLink Host PC in Camera Setup mode for participant setup/calibration
        EyelinkDoTrackerSetup(el);
    
%% 3) Qualidades dos Gabores
    % (a) Tamanho dos gabores, distância entre eles, em px, e cor do fundo
        gaborSize_px = dva2pix(screenDist,monitorW_mm/10,screenRes.width,gaborSize_dva);
        minDist_px   = dva2pix(screenDist,monitorW_mm/10, screenRes.width,minDist_dva);
        gaborBackgroundOffset = [0 0 0 0];

        minFixDist1 = gaborSize_px*fixDistFactor1;
    % (b) Cria a textura em si
        gaborTex = CreateProceduralSineGrating(window, gaborSize_px, gaborSize_px, gaborBackgroundOffset,[], 0.5);
    % (c) Propriedades do Gabor
        gaborPPD = gaborSize_px/gaborSize_dva;
        gaborFreq = gaborFreq_cpd/gaborPPD;    % Frequência, ciclos por pixel (cpp)
        gaborSigma = gaborSize_px / 6;                           % Diz respeito à nitidez da borda
                                                                 % gaborSize_px*1.5/6 mostra toda a figura, pois é maior que a diagonal 
        gaborAmplitude = .5;                                     % Como o fundo é cinza, amplitude .5 implica que vai do 0 ao 1
        gaborContrast  = 1;                                      % O contraste (entre 0 e 1) modula a amplitude
        gaborAlpha     = .5;
        phase = 0;
        gaborProps0 = [phase, gaborFreq, gaborAlpha*gaborAmplitude*gaborContrast, 0];
        gaborProps = repmat(gaborProps0', 1, nStims);
    % (d) Matrizes com orientação e centros dos estímulos. Sem os alvos, há uma 
    %     quantidade aleatória de estímulos diagonais e cardinais
        stimCenters = zeros(2, nStims, nTrials, nBlocks);
        orientation = zeros(nStims, nTrials, nBlocks);
        for b=1:nBlocks
            auxAllOri = allOri;
            auxAllOri(auxAllOri ==targetOri(b)) = [];
            orientation(:,:,b) = allOri(randi(length(auxAllOri), nStims, nTrials));
        end
        exampleGaborFactor   = 3;
        exampleGaborFreq_cpd = gaborFreq_cpd/exampleGaborFactor;
        exampleGaborSize_dva = gaborSize_dva*exampleGaborFactor;
        exampleGaborSize_px  = dva2pix(screenDist,monitorW_mm/10,screenRes.width,exampleGaborSize_dva);
        exampleGaborPPD = exampleGaborSize_px/exampleGaborSize_dva;
        exampleGaborFreq = exampleGaborFreq_cpd/exampleGaborPPD;
        exampleGaborProps = gaborProps0; exampleGaborProps(2) = exampleGaborFreq;
        exampleGaborTex = CreateProceduralSineGrating(window, exampleGaborSize_px, exampleGaborSize_px, gaborBackgroundOffset,[], 0.5);

    % (e) Cria o filtro de orientação
        [oriFilter, OFsize] = MakeOriFilter(gaborSize_px, aSigma, rSigma2);
    
%% 4) Qualidades das cruzes de fixação
        crossSize_dva = 1;
        crossSize_px = dva2pix(screenDist,monitorW_mm/10,screenRes.width,crossSize_dva);
        lineWidth_px = 4;
        fixCenters  = zeros(2, nTrials, nBlocks);
    
%% 5) Qualidades dos ruídos
    % (a) Define parâmetros do filtro e o filtro em si
        noiseCutFreq = noiseCutFreq_cpd*gaborSize_dva/gaborSize_px;
        
        noiseSigma = 1/(2*pi*noiseCutFreq); hSize = 2*ceil(noiseSigma*3)+1;
        kernel = fspecial('gaussian', hSize, noiseSigma);
    
        blurshader = EXPCreateStatic2DConvolutionShader(kernel, 1, 1, [], 1);

%         c = sqrt(log(2)) / (2 * pi);
%         noiseHiSigma = c/noiseHiCutFreq;
%         noiseLoSigma = c/noiseLoCutFreq; 
        % Usa difference of gaussians (DoG) como filtra passa-banda
%         hSize = 2*ceil(max(noiseHiSigma, noiseLoSigma)*3)+1;
%         loKernel = fspecial('gaussian', hSize, noiseLoSigma);
%         hiKernel = fspecial('gaussian', hSize, noiseHiSigma);
%         kernel = hiKernel - loKernel;    
%         blurshader = EXPCreateStatic2DConvolutionShader(kernel, 1, 1, [], 1);
    
    % (b) Define a transparência e contraste do ruído
        noiseAmplitude = .5;
        noiseAlpha = 1-gaborAlpha;% O alpha do Gabor e o alpha do ruído dizem 
                                  % respeito à contribuição deles para o resultado.
                                  % Sendo 0.5 cada, ambos contribuem igualmente.

        noiseContrast = 1;      % O contraste diz respeito ao range do ruído e
                                % do Gabor. Se for 1, o ruído assume valores
                                % do preto ao branco; se menor, o range
                                % continua simétrico ao redor de 0.5, porém
                                % mais distante dos extremos
    
    % (c) Define a matriz que contém o ruído e o centro dos retângulos a serem
    %     plotados
        noiseMatrix = zeros(gaborSize_px, nStims*gaborSize_px);
        oriPinkMatrix = zeros(gaborSize_px, nStims*gaborSize_px);
    
        noiseCenters = [0:gaborSize_px:(nStims-1)*gaborSize_px; zeros(1, nStims)] + gaborSize_px/2;
    
    % (d) Centros de cada patch de ruído
        baseRect = [0 0 gaborSize_px gaborSize_px];
        srcRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', noiseCenters(1,:), noiseCenters(2,:));
    
%% 6) Abertura Gaussiana
        blobColorOffset = [grey grey grey 0];
        [blobTex, ~] = CreateProceduralGaussBlob(window, gaborSize_px, gaborSize_px, blobColorOffset, 1, 1);
        blobContrast = 1;
        blobSigma = gaborSigma;
        blobAspect= 1;
        blobProps = [blobContrast, blobSigma, blobAspect, 0]';
    
%% 7) Localização dos estímulos (pré-computada para todos trials)
        rng('shuffle');
        for b=1:nBlocks
        for i=1:nTrials
    % (a) Cria os centros dos retângulos que contêm os Gabores e a fixação
            [currFixCenter, currStimCenter] = getStimLocations(ellipseProps, nStims, minDist_px, gridShape, randomize);
            fixCenters(:, i, b)     = currFixCenter;
            stimCenters(:, :, i, b) = currStimCenter;
    % (b) Atribui orientação não nula para apenas A(i) dos Gabores
            orientation(randperm(nStims, As(b, i)), i, b) = targetOri(b);
        end
        end
    % c) Distância mínima para considerar fixação em alvo pós-modificação
    minFixDist3 = gaborSize_px*fixDistFactor3;


%% 9) Desenha os Gabores em cada trial
    KbName('UnifyKeyNames');
    leftKey = KbName('LeftArrow');
    rightKey = KbName('RightArrow');
    spaceKey = KbName('space');
    
    try
        for b = 1:nBlocks
        for i = 1:nTrials
        % (a) Cria os retângulos de destino
            dstRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', stimCenters(1,:, i, b), stimCenters(2,:, i, b));
        
        % i. Deixa de registrar até StartRecording
            Eyelink('SetOfflineMode');
            WaitSecs(.1);

        % ii. Faz drift correction
            Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
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
                % Eyelink('Command', 'draw_circle %d %d %d %d 15', hostRect(1), hostRect(2), hostRect(3), hostRect(4));
                theta = orientation(j, i, b);
                lineLen = radius * 0.8;
                dx = lineLen * sind(theta);
                dy = lineLen * cosd(theta);
                Eyelink('Command', 'draw_line %d %d %d %d 12', round(cx-dx), round(cy-dy), round(cx+dx), round(cy+dy));
            end
        
        % (a) Desenha a cruz de fixação
            xFix = [-crossSize_px/2 crossSize_px/2 0 0];
            yFix = [0 0 -crossSize_px/2 crossSize_px/2];
            fixCoords = [xFix; yFix];
        
            Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawLines', window, fixCoords, lineWidth_px, white, fixCenters(:, i, b)', 2);
        
        % (b) Cria as texturas de ruído, sem desenhá-las
            auxNoiseMatrix = rescale(pinkNoise(gaborSize_px, nStims*gaborSize_px), -noiseAmplitude, noiseAmplitude);
            noiseMatrix(:,:) = noiseAlpha*(noiseContrast*auxNoiseMatrix)+grey;
            noiseTex = Screen('MakeTexture', window, noiseMatrix, [], [], [], [], blurshader);

        % (c) Cria as texturas dos ruídos com orientação, sem desenhá-las
            for j=1:nStims
                colRange = ((j-1)*gaborSize_px+1):(j*gaborSize_px);
                oriPinkMatrix(:,colRange) = ApplyOriFilter(oriFilter, OFsize, auxNoiseMatrix(:,colRange));
            end

            oriPinkTex = Screen('MakeTexture', window, oriPinkMatrix);

            FPonset = Screen('Flip', window);
        
        % iv. Inicia o registro da sessão
            Eyelink('StartRecording');
            Eyelink('Command', 'record_status_message "TRIAL %d/%d"', i, nTrials);
            Eyelink('Message', sprintf('FP_onset_%d', FPonset));
        
        % v. Não avança de tela até que os olhos estejam na cruz de fixação
            fixCenter = fixCenters(:,i, b);
            while true
                damn = Eyelink('CheckRecording');
                if(damn ~= 0), break; end

                if Eyelink('NewFloatSampleAvailable') > 0
                    evt = Eyelink('NewestFloatSample');
                    x_gaze = evt.gx(Eye);
                    y_gaze = evt.gy(Eye);
                    % Se os olhos estiverem perto da cruz por tempo
                    % suficiente, prossegue
                    if vecnorm([x_gaze; y_gaze] - fixCenter) <= minFixDist1
                        if (GetSecs - FPonset) >= minFixTime1, break; end
                    % Se estiverem distantes, reinicia a contagem
                    elseif vecnorm([x_gaze; y_gaze] - fixCenter) > minFixDist1
                        FPonset = GetSecs;
                    end
                end
                [~,~,keyCode] = KbCheck;
                if keyCode(KbName('space')), break; end
                WaitSecs(.01)
            end
            
            vbl = Screen('Flip', window);
            
        % (d) Desenha os ruídos, puramente opacos
            Screen('BlendFunction', window, GL_ONE, GL_ZERO);
            Screen('DrawTextures', window, noiseTex, srcRects, dstRects, orientation(:, i, b), [], [], [], []);
        
        % (e) Desenha os gratings, somando ambos os sinais
            Screen('BlendFunction', window, GL_ONE, GL_ONE);
            Screen('DrawTextures', window, gaborTex, [], dstRects, orientation(:, i, b), [], [], [], [], [], gaborProps);
        
        % (f) Desenha a abertura gaussiana
            Screen('BlendFunction', window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTextures', window, blobTex, [], dstRects, orientation(:, i, b), [], [], [0 0 0 1]', [], [], blobProps);                                                             
            Screen('Close', noiseTex);
        
            vbl = Screen('Flip', window, vbl + 0.5 * ifi);

        % vi. Registra o momento de início do trial
            trialOnset = vbl;
            Eyelink('Message', sprintf('trial_onset_%1d', trialOnset));
            Eyelink('Command', 'record_status_message "TRIAL %d', i);

        % vii. O sujeito deve visitar exatamente Ns(i) estímulos
        %      antes de ocorrerem as modificações pré-sacádicas
            counter = 0;                % Número de estímulos visitados
            flag    = zeros(1, nStims); % Quantas vezes cada estímulo foi visitado

            
            while counter < Ns(i, b)
                damn = Eyelink('CheckRecording');
                if(damn ~= 0), break; end

                if Eyelink('NewFloatSampleAvailable') > 0
                    evt = Eyelink('NewestFloatSample');
                    x_gaze = evt.gx(Eye);
                    y_gaze = evt.gy(Eye);
                    
                    % Se os olhos estiverem suficientemente próximos de
                    %  um alvo não antes visto, incrementa o contador
                    isCurrTarget = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, i, b)) <= minFixDist1;
                    if any(isCurrTarget)
                        if (GetSecs - trialOnset) >= minFixTime2
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
            end
            currIdx = find(isCurrTarget);
            seenIdx = find(flag ~= 0); notSeenIdx = find(flag == 0);


        % (g) Desenha ruído com orientação
            blinkIdx = setdiff(1:nStims, currIdx);
            Screen('BlendFunction', window, GL_ONE, GL_ZERO);
            Screen('DrawTextures', window, oriPinkTex, srcRects(:,blinkIdx), dstRects(:,blinkIdx), orientation(blinkIdx, i, b), [], [], [], []);
            Screen('BlendFunction', window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTextures', window, blobTex, [], dstRects(:,blinkIdx), orientation(blinkIdx, i, b), [], [], [0 0 0 1]', [], [], blobProps);
            Screen('Close', oriPinkTex);
    
        % viii. Registra os tempos de início e fim da 3a fase
            updateStimOnset = Screen('Flip', window);
            Eyelink('Message', sprintf('update_onset_%1d', updateStimOnset));
            updateStimOffset = Screen('Flip', window, updateStimOnset + pinkNoiseDur);
            Eyelink('Message', sprintf('update_offset_%1d', updateStimOffset));

        % ix. Verifica se há alguma fixação de duração mínima em estímulo 
        %     numa janela pós-modificação
            fixOnset = updateStimOffset; currIdx = [];
            while GetSecs - updateStimOffset < postModDur
                damn = Eyelink('CheckRecording');
                if(damn ~= 0), break; end

                if Eyelink('NewFloatSampleAvailable') > 0
                    evt = Eyelink('NewestFloatSample');
                    x_gaze = evt.gx(Eye);
                    y_gaze = evt.gy(Eye);
                    
                    isCurrTarget = vecnorm([x_gaze; y_gaze] - stimCenters(:, :, i, b)) <= minFixDist3;
                    if any(isCurrTarget)
                        if (GetSecs - fixOnset) >= minFixTime3
%                             flag = flag + isCurrTarget;
                            currIdx = find(isCurrTarget);
                            disp(['Trial ' num2str(i) ': Visitou o alvo ' num2str(currIdx) ' pós-modificação'])
                            break;
                        end
                    else
                        fixOnset = GetSecs;
                    end
                end
            end
            if ~isempty(currIdx) && ismember(currIdx, notSeenIdx)
                notSeenIdx(notSeenIdx == currIdx) = [];
            end
        
        % x. Desenha placeholders para os estímulos aleatoriamente, obedecendo
        %    a identificação deles como visitados, atual e não visitados
            allTargets = nan(1,nStims); allColors2 = allColors;
            Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            seenAux = datasample(seenIdx, min(length(seenIdx), nStimsToReport(1, i, b)), 'Replace', false);
            currAux = []; if nStimsToReport(2, i, b) == 1, currAux = currIdx; end
            notSeenAux = datasample(notSeenIdx, min(length(notSeenIdx), nStimsToReport(3, i, b)), 'Replace', false);
            orderToReportStimsCell = {seenAux, currAux, notSeenAux};
            orderToReportStims = [orderToReportStimsCell{orderToReportSets(1, i, b)} orderToReportStimsCell{orderToReportSets(2, i, b)} orderToReportStimsCell{orderToReportSets(3, i, b)}];
            
            Eyelink('Message', sprintf('answer_onset_%1d', GetSecs));
            for j=1:length(orderToReportStims)
                KbReleaseWait;
                rectColors = allColors2; rectColors(:, orderToReportStims(j)) = orange;
            
                rectPW = allPW; rectPW(:, orderToReportStims(j)) = pW2;
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
                            allColors2(:, orderToReportStims(j)) = midGrey;
                        end
                    end
            
                    allTargets(orderToReportStims(j)) = currTarget;
            
                    answerOffset = foragingFlip(window, stimCenters(:, :, i, b), dstRects, orderToReportStims, gaborSize_px, rectColors, allTargets, targetOri(b), rectPW);
            
                end
            end

            Eyelink('Message', sprintf('answer_offset_%1d', answerOffset));
            foragingFlip(window, stimCenters(:, :, i, b), dstRects, orderToReportStims, gaborSize_px, allColors, allTargets, targetOri(b), allPW);
            WaitSecs(0.5);

        % xi. Interrompe o registro, pois ou a tela será atualizada ou
        %     acabaram os trials
            Eyelink('Message', sprintf('trial_offset_%1d', GetSecs));
            Screen('Flip', window);
            WaitSecs(0.1);
            Eyelink('SetOfflineMode');
            Eyelink('StopRecording');
        end
        end
    catch
        psychrethrow(psychlasterror);
        cleanup;
    end
    Eyelink('CloseFile');
    cleanup;
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

function cleanup
    sca;
    Eyelink('Shutdown'); % Close EyeLink connection
    ListenChar(0); % Restore keyboard output to Matlab
    close all;
end

function pmf = robust_beta_pmf(m, desired_mode, concentration_type, value)
    % concentration_type: 'variance', 'sum', or 'peakness'
    % value: corresponding parameter value
    
    p = (desired_mode - 1) / (m - 1);
    
    switch concentration_type
        case 'variance'
            % Direct variance control (0 to 0.25)
            target_var = min(max(value, 0.001), 0.2);
            sol = fsolve(@(x) beta_equations(x, p, target_var), [3, 3], optimset('Display','off'));
            alpha = max(sol(1), 1.01);
            beta = max(sol(2), 1.01);
            
        case 'sum'
            % Control α+β (higher = more peaked)
            k = max(value, 2.1);
            alpha = p * (k - 2) + 1;
            beta = k - alpha;
            
        case 'peakness'
            % Simple intuitive control (1-10 scale)
            k = 2 + value * 8;  % Map to [2, 10]
            alpha = p * (k - 2) + 1;
            beta = k - alpha;
    end
    
    % Generate distribution
    x_continuous = (0.5:(m-0.5)) / m;
    densities = betapdf(x_continuous, alpha, beta);
    pmf = densities / sum(densities);
end

function currTime = foragingFlip(win, centers, dstCoord, idxForOvals, s, colors, aT, tOri, pW)
    Screen('FrameOval', win, colors(:, idxForOvals), dstCoord(:, idxForOvals), pW(:, idxForOvals));

    targets = (aT == true); nonTargets = (aT == false);

    drawInclinedLines(win, centers(:,targets), s*.8, tOri, colors(:,targets), pW(:, targets));
    drawDots(win, centers(:,nonTargets), s/5, colors(:,nonTargets));

    currTime = Screen('Flip', win);
    WaitSecs(0.001);
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

function drawDots(windowPtr, centers, len, color)
    if ~isempty(centers)
        Screen('DrawDots', windowPtr, centers, len, color, [], 2);
    end
end


function runDemo(win)

    DrawFormattedText(win, ...
        'Demo running...\n\nPress SPACE to return to menu', ...
        'center', 'center', WhiteIndex(win));
    Screen('Flip', win);
    
    while true
        [~,~,keyCode] = KbCheck;
        if keyCode(KbName('space'))
            KbReleaseWait;
            break;
        end
    end

end

function runTraining(win)

DrawFormattedText(win, ...
    'Training mode\n\nPress SPACE to return', ...
    'center', 'center', WhiteIndex(win));
Screen('Flip', win);

KbStrokeWait;
end
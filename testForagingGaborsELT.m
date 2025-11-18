function testForagingGaborsELT(nTrials, nStims, nTargets)
% Versão com Gabores + Eyelink
% Parâmetros dos estímulos que podem ser ajustados:
% gaborAlpha: o quanto o gabor contribui em relação ao ruído
% gaborContrast: modula a amplitude do Gabor
% noiseContrast: modula a amplitude do ruído

        if nargin < 1, nTrials = 10; end
        if nargin < 2, nStims   = 8; end
        if nargin < 3, nTargets = 1; end
    
    
%% 1) Parâmetros das condições experimentais
    % (a) Distância da tela, em cm
        screenDist       = 57; % PC: 57; tv: 200

    % (b) Tamanhos e frequências, em dva e cpd
        gaborSize_dva    = 1;
        gaborFreq_cpd    = 11; % usa 11
        noiseCutFreq_cpd = 5; % em cpd

    % (c) Disposição básica dos estímulos
        minDist_dva = max(3, gaborSize_dva);% OBS: ajustar com ellipseToScreenRatio 
        ellipseToScreenRatio = [3/4 3/4];   % Proporção da elipse em relação ao
                                              % tamanho da tela em cada dimensão
        gridShape   = [2 4];                % #linhas x #colunas
        randomize = false;
    % (d) Tempo mínimo de fixação na cruz
        minFixTime = .5; % em s
        fixDistFactor = 1.2;
    
    
%% 2) Inicializa PTB e EyelinkToolBox
    % i. Caixa de diálogo
        prompt = {'Subject code', 'Session number', 'Dominant eye (L or R)'};
        dlg_title = 'Experimental Setup';
        def = {'00', '01', 'R'};
        options.Resize = 'on';
        answer = inputdlg(prompt, dlg_title, 1, def, options);
        domEye = answer{3};
        if  isempty(answer), fprintf('Session cancelled by user\n'); cleanup; return; end
        edfFile = [answer{1} '_' answer{2}];

        % Interrompe se nome muito longo
        if length(edfFile) > 8, fprintf('Filename needs to be no more than 8 characters long (letters, numbers and underscores only)\n'); cleanup; return; end
    
    
    % (a) Configurações antes de abrir a tela
        FlushEvents;
        PsychDefaultSetup(2);
        Screen('Preference', 'SyncTestSettings', 0.01, 50, 0.25);
        Screen('Preference', 'SuppressAllWarnings', 1);
        Screen('Preference', 'Verbosity', 0);
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    

    % (b) Escolhe a tela em que haverá o desenho
        screenNumber = max(Screen('Screens'));
        white = WhiteIndex(screenNumber);
        grey = white / 2;%startP startP startP+dimX startP+dimY
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

    % i. Inicializa Eyelink
        dummymode = 0;
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

    % i. Atualiza o tamanho mínimo da fixação
        minFixDist = gaborSize_px*fixDistFactor;
    % (b) Cria a textura em si
        gaborTex = CreateProceduralSineGrating(window, gaborSize_px, gaborSize_px, gaborBackgroundOffset,[], 0.5);
    % (c) Propriedades do Gabor
        gaborFreq = gaborFreq_cpd*gaborSize_dva/gaborSize_px;    % Frequência, ciclos por pixel (cpp)
        gaborSigma = gaborSize_px / 6;                           % Diz respeito à nitidez da borda
                                                                 % gaborSize_px*1.5/6 mostra toda a figura, pois é maior que a diagonal 
        gaborAmplitude = .5;                                     % Como o fundo é cinza, amplitude .5 implica que vai do 0 ao 1
        gaborContrast  = 1;                                      % O contraste (entre 0 e 1) modula a amplitude
        gaborAlpha     = .5;
        phase = 0;
        gaborProps = [phase, gaborFreq, gaborAlpha*gaborAmplitude*gaborContrast, 0];
        gaborProps = repmat(gaborProps', 1, nStims);
    % (d) Matrizes com orientação e centros dos estímulos
        stimCenters = zeros(2, nStims, nTrials);
        orientation = 45*ones(nTrials, nStims);
    
%% 4) Qualidades das cruzes de fixação
        lineWidth_px = 4;
        fixCenters  = zeros(2, nTrials);
    
%% 5) Qualidades dos ruídos
    % (a) Define parâmetros do filtro e o filtro em si
        noiseCutFreq = noiseCutFreq_cpd*gaborSize_dva/gaborSize_px;
        
        noiseSigma = 1/(2*pi*noiseCutFreq); hSize = 2*ceil(noiseSigma*3)+1;
        kernel = fspecial('gaussian', hSize, noiseSigma);
    
        blurshader = EXPCreateStatic2DConvolutionShader(kernel, 1, 1, [], 1);
    
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
        for i=1:nTrials
    % (a) Cria os centros dos retângulos que contêm os Gabores e a fixação
            [currFixCenter, currStimCenter] = getStimLocations(ellipseProps, nStims, minDist_px, gridShape, randomize);
            fixCenters(:,i) = currFixCenter;
            stimCenters(:,:, i) = currStimCenter;
    % (b) Atribui orientação não nula para apenas nTargets dos Gabores
            orientation(i, randperm(nStims, nTargets)) = 90;
        end
        
    
%% 8) Desenha os Gabores em cada trial
    KbName('UnifyKeyNames');
    
    try
        for i = 1:nTrials
        % (a) Cria os retângulos de destino
            dstRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', stimCenters(1,:, i), stimCenters(2,:, i));
        
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
                theta = orientation(i, j);
                lineLen = radius * 0.8;
                dx = lineLen * sind(theta);
                dy = lineLen * cosd(theta);
                Eyelink('Command', 'draw_line %d %d %d %d 12', round(cx-dx), round(cy-dy), round(cx+dx), round(cy+dy));
            end
        
        % (a) Desenha a cruz de fixação
            xFix = [-gaborSize_px/2 gaborSize_px/2 0 0];
            yFix = [0 0 -gaborSize_px/2 gaborSize_px/2];
            fixCoords = [xFix; yFix];
        
            Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawLines', window, fixCoords, lineWidth_px, white, fixCenters(:,i)', 2);
        
        % (b) Cria as texturas de ruído, sem desenhá-las
            noiseMatrix(:,:) = noiseAlpha*(noiseContrast*rescale(pinkNoise(gaborSize_px, nStims*gaborSize_px), -noiseAmplitude, noiseAmplitude))+grey;
            noiseTex = Screen('MakeTexture', window, noiseMatrix, [], [], [], [], blurshader); 
            
            FPonset = Screen('Flip', window);
        
        % iv. Inicia o registro da sessão
            Eyelink('StartRecording');
            Eyelink('Command', 'record_status_message "TRIAL %d/%d"', i, nTrials);
            Eyelink('Message', sprintf('FP_onset_%d', FPonset));
        
        % v. Não avança de tela até que os olhos estejam na cruz de fixação
            fixCenter = fixCenters(:,i);
            while true
                damn = Eyelink('CheckRecording');
                if(damn ~= 0), break; end

                if Eyelink('NewFloatSampleAvailable') > 0
                    evt = Eyelink('NewestFloatSample');
                    x_gaze = evt.gx(Eye);
                    y_gaze = evt.gy(Eye);
                    % Se os olhos estiverem perto da cruz por tempo
                    % suficiente, prossegue
                    % vecnorm([x_gaze; y_gaze] - fixCenter)
                    % minFixDist
                    if vecnorm([x_gaze; y_gaze] - fixCenter) <= minFixDist
                        if (GetSecs - FPonset) >= minFixTime, break; end
                    % Se estiverem distantes, reinicia a contagem
                    elseif vecnorm([x_gaze; y_gaze] - fixCenter) > minFixDist
                        FPonset = GetSecs;
                    end
                end
                % [~,~,keyCode] = KbCheck;
                % if keyCode(KbName('space')), break; end
            end
            
            vbl = Screen('Flip', window);
            
        % (d) Desenha os ruídos, puramente opacos
            Screen('BlendFunction', window, GL_ONE, GL_ZERO);
            Screen('DrawTextures', window, noiseTex, srcRects, dstRects, orientation(i,:), [], [], [], []);
        
        % (e) Desenha os gratings, somando ambos os sinais
            Screen('BlendFunction', window, GL_ONE, GL_ONE);
            Screen('DrawTextures', window, gaborTex, [], dstRects, orientation(i,:), [], [], [], [], [], gaborProps);
        
        % (f) Desenha a abertura gaussiana
            Screen('BlendFunction', window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTextures', window, blobTex, [], dstRects, orientation(i,:), [], [], [0 0 0 1]', [], [], blobProps);                                                             
            Screen('Close', noiseTex);
        
            vbl = Screen('Flip', window, vbl + 0.5 * ifi);

        % vi. Registra o momento de início do trial
            trialOnset = vbl;
            Eyelink('Message', sprintf('trial_onset_%1d', trialOnset));
            Eyelink('Command', 'record_status_message "TRIAL %d', i);

        % vii. Quando a barra de espaço for apertada, o trial é encerrado
            while true
                [~,trialOffset,keyCode] = KbCheck;
                if keyCode(KbName('space')), break; end
            end
            Eyelink('Message', sprintf('trial_offset_%1d', trialOffset));
            Screen('Flip', window);
         
        % vii. Interrompe o registro, pois ou a tela será atualizada ou
             % acabaram os trials
            WaitSecs(0.1); 
%             Eyelink('SetOfflineMode');
%             Eyelink('StopRecording');

            disp(['Termina o ' num2str(i) 'o trial (após um Flip cinza)'])
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
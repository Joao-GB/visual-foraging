function testForagingGabors(nTrials, nStims, nTargets)
% Por enquanto, a função apenas produz nStims Gabores, apresentados ao 
% longo de nTrials trials, e, dentro de cada trial, há exatamente nTargets
% alvos.

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
    
    
%% 2) Inicializa PTB
    % (a) Configurações antes de abrir a tela
        PsychDefaultSetup(2);
        Screen('Preference', 'SkipSyncTests', 2);
        Screen('Preference','SuppressAllWarnings',0);
        Screen('Preference','Verbosity',0);
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    
    % (b) Escolhe a tela em que haverá o desenho
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
    
    
%% 3) Qualidades dos Gabores
    % (a) Tamanho dos gabores, distância entre eles, em px, e cor do fundo
        gaborSize_px = dva2pix(screenDist,monitorW_mm/10,screenRes.width,gaborSize_dva);
        minDist_px   = dva2pix(screenDist,monitorW_mm/10, screenRes.width,minDist_dva);
        gaborBackgroundOffset = [0 0 0 0];
    % (b) Cria a textura em si
          gaborTex = CreateProceduralSineGrating(window, gaborSize_px, gaborSize_px, gaborBackgroundOffset,[], 0.5);
    % (c) Propriedades do Gabor
        gaborFreq = gaborFreq_cpd*gaborSize_dva/gaborSize_px;    % Frequência, ciclos por pixel (cpp)
        gaborSigma = gaborSize_px / 6;                           % Diz respeito à nitidez da borda
                                                                 % gaborSize_px*1.5/6 mostra toda a figura, pois é maior que a diagonal 
        gaborAmplitude = .5;                                     % Como o fundo é cinza, amplitude .5 implica que vai do 0 ao 1
        gaborContrast  = .5;                                      % O contraste (entre 0 e 1) modula a amplitude
        gaborAlpha     = .5;
        phase = 0;
        gaborProps = [phase, gaborFreq, gaborAlpha*gaborAmplitude*gaborContrast, 0];
        gaborProps = repmat(gaborProps', 1, nStims);
    % (d) Matrizes com orientação e centros dos estímulos
        stimCenters = zeros(2, nStims, nTrials);
        orientation = zeros(nTrials, nStims);
    
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

        noiseContrast = .5;      % O contraste diz respeito ao range do ruído e
                                % do Gabor. Se for 1, o ruído assume valores
                                % do preto ao branco; se menor, o range
                                % continua simétrico ao redor de 0.5, porém
                                % mais distante dos extremos
    
    % (c) Define a matriz que contém o ruído e o centro dos retêngulos a serem
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
    vbl = Screen('Flip', window);
    
    for i = 1:nTrials
    % (a) Desenha a cruz de fixação
        xFix = [-gaborSize_px/2 gaborSize_px/2 0 0];
        yFix = [0 0 -gaborSize_px/2 gaborSize_px/2];
        fixCoords = [xFix; yFix];
    
        Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('DrawLines', window, fixCoords, lineWidth_px, white, fixCenters(:,i)', 2);
    
    % (b) Cria as texturas de ruído, sem desenhá-las
        noiseMatrix(:,:) = noiseAlpha*(noiseContrast*rescale(pinkNoise(gaborSize_px, nStims*gaborSize_px), -noiseAmplitude, noiseAmplitude))+grey;
        noiseTex = Screen('MakeTexture', window, noiseMatrix, [], [], [], [], blurshader); 
    
        Screen('Flip', window, vbl + 0.5 * ifi);
    
        while true
            [~,~,keyCode] = KbCheck;
            if keyCode(KbName('space')), break; end
        end
    
    % (c) Cria os retângulos de destino
        dstRects = CenterRectOnPointd(repmat(baseRect, [nStims,1])', stimCenters(1,:, i), stimCenters(2,:, i));
        size(dstRects)
        
        vbl = Screen('Flip', window);
        
    % (d) Desenha os ruídos, puramente opacos
        Screen('BlendFunction', window, GL_ONE, GL_ZERO);
        Screen('DrawTextures', window, noiseTex, srcRects, dstRects, [], [], [], [], []);
    
        Screen('BlendFunction', window, GL_ONE, GL_ONE);
    % (e) Desenha os gratings, somando ambos os sinais
        Screen('DrawTextures', window, gaborTex, [], dstRects, orientation(i,:), [], [], [], [], [], gaborProps);
        
    
    % (f) Desenha a abertura gaussiana
        Screen('BlendFunction', window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
        Screen('DrawTextures', window, blobTex, [], dstRects, [], [], [], [0 0 0 1]', [], [], blobProps);                                                             
        Screen('Close', noiseTex);
    
        vbl = Screen('Flip', window, vbl + 0.5 * ifi);
        KbStrokeWait;
    
        
    end
        sca;
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

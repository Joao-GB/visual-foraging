function params = foragingParams
%% Propriedades do setup experimental
    % (a) Distância entre sujeito e tela, em cm
    params.screenDist    = 57;

%% Propriedades de aparência dos estímulos
    % (b) Parâmetros das cruzes de fixação
    params.crossSize_dva = 1;
    params.lineWidth_px = 4;

    % (c) Tamanhos, frequência e outras propriedades dos gabores
    params.gaborSize_dva  = 1.5;
    params.gaborFreq_cpd  = 15; % Rucci et al. (2007) usam 11, em cpd
    params.gaborAmplitude = .5;
    params.gaborContrast  = 1;
    params.gaborAlpha     = .5;
    params.gaborPhase     = 0;

    % (d) Parâmetros do ruído
    params.noiseAmplitude = .5;
    params.noiseAlpha     = 1 - params.gaborAlpha;
    params.noiseContrast  = 1;

    % (e) Frequências de corte para filtro a ser aplicado no ruído, em cpd
    params.noiseLoCutFreq_cpd = 1/params.gaborSize_dva;
    params.noiseHiCutFreq_cpd = max(5,2*params.noiseLoCutFreq_cpd);

    % (f) Parâmetros do ruído rosa com orientação
    %% Arrumar parâmetros
    params.aSigma = 20;
    params.rSigma2 = .76;

    % (g) Parâmetros dos blobs
    params.blobContrast = 1;
    params.blobAspect   = 1;

    % (h) Regras de disposição dos estímulos
        % OBS: Ajustar minDist_dva junto com ellipseToScreenRatio
    params.minDist_dva = max(3, params.gaborSize_dva); 
    params.ellipseToScreenRatio = [3/4 3/4];
    params.gridShape   = [2 4];
    params.randomize = false;

    % (i) Orientações possíveis dos estímulos, em graus
    params.allOri = [0 45 90 135]; %[0 90];%
    params.allOriMap  = containers.Map(params.allOri,1:length(params.allOri));
    params.allOriName = {'vertical', 'diagonal crescente', 'horizontal', 'diagonal decrescente'};%{'vertical', 'horizontal'};%

%% Propriedades indispensáveis da tarefa: tempos de fixação, durações dos estímulos...
    % (j) Quantidades e proporções de estímulos a serem reportados ao 
    %     longo dos trials
    params.minToReport = 2;
    params.maxToReport = 3;
    params.propTrialsPSA = .8;                 % Proporção de trials em que se deve reportar estímulo fixado pós-modificação

    % (k) Duração do estímulo (ruído rosa), em s
    params.pinkNoiseDur = .15;

    % (l) Parâmetros temporais de fixação
    params.minFixTime1 = .5;                    % Tempo mínimo de fixação na cruz inicial, em s
    params.minFixTime2 = .1;                   % Tempo mínimo de fixação nos estímulos, em s
    params.medFixTime2 = .3;                    % Tempo médio de fixação nos estímulos, em s, a ser usado apenas para a fila inicial de fixações
    params.minFixTime3 = .1;                   % Tempo mínimo de fixação na região do alvo pós-modificação, em s
    params.postModDur  = .75;                    % Janela temporal (após offset dos estímulos) durante a qual se 
                                                    % espera a fixação com duração mínima minFixTime3, em s
    params.blobPMDur   = .05;
    % (m) Tamanho da fila de tempos de fixação
    params.fixTimeQueueSize = 20;

    % (n) Fatores multiplicativos para se considerar ou não fixação
    params.fixDistFactor1 = 1;    % Fator de tolerância para distância entre fixação e alvos pré-modificação
    params.fixDistFactor3 = 1.3;    % Fator de tolerância para distância entre fixação e alvos pós-modificação
                                        % (maior tolerância a erro, já que o estímulo é removido)


%% -----
%% Propriedades secundárias da tarefa: treino e cursor, tempos limites, retentativas...
    % (a) Número de trials no treino
    params.nTrialsTrain  = 5;

    % (b) Tempo máximo durante o qual é exibida a cruz antes de mensagem de
    % 'Tempo Esgotado' na fase 1
    params.maxCrossDur = 6;

    % (c) Fatores de duração de cada trial com base nas durações de fixações e sacadas
    params.revisitFactor = 1/2;         % Proporção de tolerância de revisita a
                                            % estímulos já vistos
    params.sacFixDurRatio = 1;        % Duração de uma sacada em relação a uma fixação
    params.maxTrialDurFactor = (1+params.revisitFactor)*(1+params.sacFixDurRatio); 
                                            % A ideia é que se N é um número de estímulos, o sujeito fará no máximo
                                            % 1.5*N fixações e sacadas, já que a 1a fixação parte da cruz de fixação.
                                            % Por exemplo, se revisitFactor = 1/2 e sacFixDurRatio = 1/2, sendo m o 
                                            % tempo de cada fixação (logo m/2 de cada sacada), o tempo máximo é 
                                            % 1.5(N*m + N*(.5*m)) = 2.25*N*m
    params.cursorSacFixDurRatio = 5;    % Mesma finalidade que sacFixDurRatio,
                                            % mas para demonstração com cursor
    params.cursorMaxTrialDurFactor = (1+params.revisitFactor)*(1+params.cursorSacFixDurRatio);

    params.cursorPinkNoiseDur   = .25;
    params.cursorPostModDur     = 1;

    % (d) Parâmetros de retentativas
    params.nBufferTrials = 3;   % Número de trials adicionados para evitar
                                    % repetições consecutivas de um trial
    params.maxRetries = 3;      % Número máximo de vezes que pode se 
                                    % repetir cada trial 
    % (e) Tamanho e cor do cursor
    params.cursorRadius_px   = 4;
    params.cursorColor       = [1 1 1];

    % (f) Parâmetros de fixação com cursor
    params.cursorMinFixTime2 = .2;
    params.cursorMedFixTime2 = .5;
    params.cursorMinFixTime3 = params.cursorMinFixTime2/2;

    % (g) Textos de início, fim, pausa e interrupção de sessão, blocos e
    %     trials
    params.msg.suffix     = {'cur', 'trn', 'exp'};
    params.msg.on.ses{1}     = 'SESSION ONSET (%s)';
    params.msg.on.ses{2}     = 'SESSION ONSET';
    params.msg.on.blk{1}  = 'BLOCK ONSET %d/%d';
    params.msg.on.blk{2}  = 'BLOCK ONSET';
    params.msg.on.trl{1}  = 'TRIAL ONSET %d/%d';
    params.msg.on.trl{2}  = 'TRIAL ONSET';
    params.msg.on.stm{1}  = 'STIM ONSET NEW';
    params.msg.on.stm{2}  = 'STIM ONSET OLD';
    params.msg.on.P1      = 'PHASE 1 ONSET';
    params.msg.on.P2      = 'PHASE 2 ONSET';
    params.msg.on.P3      = 'PHASE 3 ONSET';
    params.msg.on.PM      = 'POST-MOD ONSET';
    params.msg.on.P4      = 'PHASE 4 ONSET';
    params.msg.on.pse     = 'PAUSE ONSET';
    

    params.msg.err.blk    = 'BLOCK ABORT TRY';
    params.msg.err.trl{1} = 'TRIAL ABORT TRY';
    params.msg.err.trl{2} = 'TRIAL ABORT IDX';
    params.msg.err.P1     = 'TRIAL TIMEOUT PHASE 1';
    params.msg.err.P2     = 'TRIAL TIMEOUT PHASE 2';

    params.msg.pse{1}     = 'PAUSE > MENU';
    params.msg.pse{2}     = 'PAUSE > RESTART';
    params.msg.pse{3}     = 'PAUSE > RECAL';
    params.msg.pse{4}     = 'PAUSE > RESUME';
    params.msg.pse{5}     = 'TRIAL > RECAL';

    params.msg.off.ses{1} = 'SESSION OFFSET (%s) COMPLETED';
    params.msg.off.ses{2} = 'SESSION OFFSET (%s) INTERRUPTED';
    params.msg.off.ses{3} = 'SESSION OFFSET';
    params.msg.off.blk{1} = 'BLOCK OFFSET %d/%d';
    params.msg.off.blk{2} = 'BLOCK OFFSET';
    params.msg.off.trl{1} = 'TRIAL OFFSET %d/%d';
    params.msg.off.trl{2} = 'TRIAL OFFSET';
    params.msg.off.stm{1} = 'STIM OFFSET';
    params.msg.off.stm{2} = 'STIM OFFSET BAD';
    params.msg.off.stm{3} = 'STIM OFFSET P3';
    params.msg.off.stm{4} = 'STIM OFFSET PM';
    params.msg.off.P3     = 'PHASE 3 OFFSET';
    params.msg.off.PM      = 'POST-MOD OFFSET';
    params.msg.off.P4     = 'PHASE 4 OFFSET';
    params.msg.off.pse    = 'PAUSE OFFSET';

%% Parâmetros de formatação: fontes, espessuras, cores, tamanhos e posições de textos e ícones...
    % (g) Fonte para os textos
    params.textFont          = 'Roboto';

    % (h) Diferentes tamanhos de fonte
    params.textSizeTitle     = 110;
    params.textSizeEnormous  = 90;
    params.textSizeHuger     = 80;
    params.textSizeHuge      = 60;
    params.textSizeBigger    = 45;
    params.textSizeBig       = 30;
    params.textSizeNormal    = 24;
    params.textSizeNormalish = 20;
    params.textSizeSmall     = 14;
    params.textSizeSmaller   = 10;
    params.textSizeTiny      = 8;
    params.textSizeTinier    = 5;

    % (i) Espessura de caneta de desenho
    params.pW1 = 3;
    params.pW2 = 5;
    params.pW3 = 10;

    % (j) Cores escolhidas a dedo
    params.orange   = [.98 .55 .06];
    params.red      = [.86 .20 .20];
    params.blue     = [.57 .64 .72];
    params.darkBlue = [.40 .50 .60];
    params.green    = [.06 .87 .33];

    % (k) Tamanhos e posições dos símbolos de Menu e Pause
    params.titleMarginFactor = .05;     % Distância entre título e limites da tela

    params.btnW   = 300;                % Tamanho dos botões clicáveis
    params.btnH   = 300;

    params.gap    = 80;                 % Distância máxima entre botões
    params.upFrac = 2/3;                % Proporção do botão reservada para ícones
    params.iconScaleFactor = 0.6;       % Proporção da porção reservada realmente 
                                            % ocupada
    params.outFolder       = 'out';     % Nome da pasta com imagens dos ícones
    params.imgFolder       = 'img';     % Nome da pasta com imagens dos ícones
    params.imgExtension    = '.png';    % Extensão das imagens

%% Durações e delays de detalhes
    params.repeatDelay = 0.5;   % Tempo que demora para considerar um segundo
                                    % aperto de tecla (se for mantida apertada)
    params.repeatRate  = 0.1;   % Tempo com que uma tecla é reconsiderada 
                                    % apertada passado o repeatDelay
    params.pauseMaxDur = 15;    % Duração máxima de uma pausa sem se seja exibida
                                    % tela relembrando alvo
    % (l) Atrasos para começar efeitos de fade in
    params.fadeInDelay1 = 0.3;
    params.fadeInDelay2 = 0.7;

    % (m) Durações de efeitos de fade in
    params.fadeInDur1   = 0.1;
    params.fadeInDur2   = 0.35;
    params.fadeInDur3   = 1.5;
end
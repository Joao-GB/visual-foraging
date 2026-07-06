function [trlProps, analysis, eyeData, evTimes] = foragingAnalysis(subj, ses)
    % Convenção:
    % vou chamar os estímulos forrageado (for), sacádico (sacc) e 
    % não-sacádico (nSacc) e identificá-los como -1, 0 e 1 com respeito ao 
    % tempo em que são vistos quanto ao início da sacada

    % Adiciona os caminhos necessários
    currFolder = fileparts(mfilename('fullpath')); parentFolder = fileparts(currFolder);
    addpath(genpath(fullfile(currFolder, 'dep')));
    addpath(genpath(fullfile(currFolder, 'plt')));
    addpath(parentFolder);
    params = foragingParams; outFolder = fullfile(params.currFolder, params.outFolder);

    % Carrega os arquivos da sessão
    [mat, edf, sesStr] = foragingLoad(subj, ses, outFolder, params);
    mat.results.trialOrder(2,:,:) = mat.results.trialOrder(2,:,:) .* reshape(1:mat.tkP.nBlocks, 1, 1, []);

    % Extrai as informações de cada trial, que já passaram por uma pré-seleção
    [trlProps, eyeData, evTimes] = foragingTrlProps(mat, edf, sesStr, subj);
    
    [allGoodTrl, ~] = getGoodTrl(trlProps, mat);
    trlProps = trlProps(allGoodTrl);

    clear currFolder edf outFolder params parentFolder sesStr ses subj;

    %% -- TAREFA PRÉ-SACÁDICA --
    % 1. Análise exploratória de estímulos e comportamento ocular na tarefa pré-sacádica
    preProbePos        = pixel_to_dva([trlProps.preProbePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    preProbePosFix     = pixel_to_dva([trlProps.preProbePosFixPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    probePos      = pixel_to_dva([trlProps.probePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    probePosFix   = pixel_to_dva([trlProps.probePosFixPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    nSaccProbePos = pixel_to_dva([trlProps.nSaccProbePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';

      % (a) Posição relativa entre pré-probe, probe e nSacc probe
    plotPSAStimPos(preProbePos, preProbePosFix, probePos, probePosFix, nSaccProbePos, mat.drP);

      % (b) Caracterização do triângulo definido por esses três pontos
    plotPSAStimTriangleProps(preProbePos, probePos, nSaccProbePos, mat.drP);

      % (c) Direção e amplitude das sacadas
    plotPSASaccProps(preProbePos, probePos, preProbePosFix, probePosFix, mat.drP);

    % 2. Gráfico principal: Efeito pré-sacádico E tabela de contingência
    PSA = plotPSAmain(trlProps, mat.drP);

    % 3. Efeito da ordem das perguntas
    order = plotPSAorder(trlProps, mat.drP);
    PSA.order = order;

    % Efeito de categorias do pré-probe e do probe
    %% Adicionar d-prime
    plotPSAcat([trlProps.preProbeCat], [trlProps.probeCat], [trlProps.probeHit], [trlProps.nSaccProbeHit], mat.drP);

    % Efeito do desempenho no forrageamento
    %% Adicionar d-prime
    plotPSAforagingPerformance(trlProps, PSA, mat.drP);
    
    % Efeito do número de vistos (no pré-s e na duração da fixação)
    %% Adicionar d-prime
    plotPSAforagingNumSeen(trlProps, mat.drP)

    % Durações relevantes: ruído rosa, tempo de fixação e latência da sacada pós ruído rosa
    plotPSAdurations(trlProps, mat.drP)

    % Efeito da duração do ruído rosa (split)
    %% Adicionar d-prime
    plotPSApinkDur(trlProps, mat)

    % Efeito da duração da fixação no desempenho
    %% Adicionar d-prime
    plotPSAfixDur(trlProps, mat)

    %% Efeito da ditância na performance (idealmente, invariante para pré-s mas variável para n-sac)
    %% Adicionar d-prime
    % Talvez ordenar e usar quantis de modo que todos bins tenham mesma quantidade de pontos, e a distância média representada por eles?

    %% Mapa de d primes em função da distância e orientação, em relação a nSacc e probe
    % Alinhar tanto com posição destino como posição final da fixação
    plotPSAdprimeMap


    %% -- TAREFA DE FORRAGEAMENTO --
    %% Efeito do número de fixações anteriores na duração da atual

    %% Efeito do número de vistos no desempenho do forragemaento

    %% Efeito da distribuição espacial ou temporal (i.e., quão longe ou há quanto tempo faz que foi visto) no desempenho do forrageamento
    % Ou quão isolado

    %% Duração da fixação em função de ser alvo ou não

    %% -- PROPRIEDADES DA TAREFA --
    %% Distância mínima entre estímulos


end
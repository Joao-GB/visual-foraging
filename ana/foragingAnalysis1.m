function [trlProps, analysis, eyeData, evTimes] = foragingAnalysis1(subj, ses, searchFolder)
    % Convenção:
    % vou chamar os estímulos forrageado (for), sacádico (sacc) e 
    % não-sacádico (nSacc) e identificá-los como -1, 0 e 1 com respeito ao 
    % tempo em que são vistos quanto ao início da sacada

    % Adiciona os caminhos necessários
    currFolder = fileparts(mfilename('fullpath')); parentFolder = fileparts(currFolder);
    addpath(genpath(fullfile(currFolder, 'dep')));
    addpath(genpath(fullfile(currFolder, 'plt')));
    addpath(parentFolder);
    params = foragingParams; 
    if nargin < 3, outFolder = fullfile(params.currFolder, params.outFolder);
    else,          outFolder = searchFolder;
    end

    % Carrega os arquivos da sessão
    [mat, edf, sesStr] = foragingLoad(subj, ses, outFolder, params);
    mat.results.trialOrder(2,:,:) = mat.results.trialOrder(2,:,:) .* reshape(1:mat.tkP.nBlocks, 1, 1, []);

    % Extrai as informações de cada trial, que já passaram por uma pré-seleção
    [trlProps, eyeData, evTimes] = foragingTrlProps(mat, edf, sesStr, subj);
    
    [allGoodTrl allTrl] = getGoodTrl(trlProps, mat, 'TRIALS EXPERIMENTO');
    auxTrlProps = trlProps;
    trlProps = trlProps(allGoodTrl);

%     if isfield(mat.tkP.t1.tkP, 'nTrainTrials')
%         [trainTrlProps, ~, ~] = foragingTrlProps(mat, edf, sesStr, subj, 2);
%         [trainAllGoodTrl, trainAllTrl] = getGoodTrl(trainTrlProps, mat, 'TRIALS TREINO');
%         trainTrlProps1 = trainTrlProps(trainAllGoodTrl);
%         trainTrlProps2 = trainTrlProps(trainAllTrl);
%     
%         trainTrlProps = trainTrlProps1;
%     end

    clear currFolder edf outFolder params parentFolder sesStr ses subj;

    skipAll = false;

    %% -- TAREFA PRÉ-SACÁDICA --
    fprintf('\n\n-- TAREFA PRÉ-SACÁDICA --\n')

    % 1. Análise exploratória de estímulos e comportamento ocular na tarefa pré-sacádica
    preProbePos        = pixel_to_dva([trlProps.preProbePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    preProbePosFix     = pixel_to_dva([trlProps.preProbePosFixPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    probePos      = pixel_to_dva([trlProps.probePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    probePosFix   = pixel_to_dva([trlProps.probePosFixPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    nSaccProbePos = pixel_to_dva([trlProps.nSaccProbePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';

    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Gráficos de caracterização espacial e temporal da tarefa pré-sacádica');
        
        auxPlot = strcmp(skipMode, 'none');
        if closePlots, close all; end

          % (a) Posição relativa entre pré-probe, probe e nSacc probe
        [rotProbePos, rotProbeFix, rotNSaccProbePos] = plotPSAStimPos(preProbePos, preProbePosFix, probePos, probePosFix, nSaccProbePos, mat.drP, auxPlot);
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            
            
            disp('> Gráficos de caracterização espacial e temporal da tarefa pré-sacádica');
        
              % (b) Caracterização do triângulo definido por esses três pontos
            plotPSAStimTriangleProps(preProbePos, probePos, nSaccProbePos, mat.drP);
        
              % (c) Direção e amplitude das sacadas
            plotPSASaccProps(preProbePos, probePos, preProbePosFix, probePosFix, mat.drP);

              % (d) Durações relevantes: ruído rosa, tempo de fixação e latência da sacada pós ruído rosa
            plotPSAdurations(trlProps, mat.drP)
        end
    end

    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Gráficos de efeito PSA principal');
        
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Gráficos de efeito PSA principal');
              % 2. Gráfico principal: Efeito pré-sacádico E tabela de contingência
            PSA = plotPSAmain(trlProps, mat.drP);
            sHit = PSA.main.sacc.idx{1,1}+PSA.main.sacc.idx{2,2};
            nHit = PSA.main.nSacc.idx{1,1}+PSA.main.nSacc.idx{2,2};
        
            %% Para o treino:
%             if isfield(mat.tkP.t1.tkP, 'nTrainTrials')
%                 trainPSA = plotPSAmain(trainTrlProps, mat.drP, 1);
%             end
        end
    end

    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Gráficos de efeito PSA secundários');
        
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Gráficos de efeito PSA secundários');

            % Efeito de categorias do pré-probe e do probe
            plotPSAcat([trlProps.preProbeCat], [trlProps.probeCat], [trlProps.probeHit], [trlProps.nSaccProbeHit], mat.drP);
            
            % Efeito da direção da sacada (tanto em ângulo como em 'faixa' de
            % amplitude ao redor da direção)
            plotSacDirEffect(rotProbeFix, rotProbePos, rotNSaccProbePos, sHit, nHit, mat.drP)
            plotSacRangeEffect(rotProbeFix, rotProbePos, rotNSaccProbePos, sHit, nHit, mat.drP)
        end
    end

    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Gráficos de efeito PSA terciários');
        
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Gráficos de efeito PSA terciários');
              % 3. Efeito da ordem das perguntas
            order = plotPSAorder(trlProps, mat.drP);
            PSA.order = order;
        end
    end

    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Gráficos de efeito PSA e forrageamento');
        
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Gráficos de efeito PSA e forrageamento');
            % Efeito do desempenho no forrageamento
            plotPSAforagingPerformance(trlProps, PSA, mat.drP);
            
            % Efeito do número de vistos (no pré-s e na duração da fixação)
            plotPSAforagingNumSeen(trlProps, mat.drP)
        end
    end

    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Gráficos de efeito PSA vs durações e distâncias');
        
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Gráficos de efeito PSA vs durações e distâncias');
            % Efeito da duração do ruído rosa (split)
            plotPSApinkDur(trlProps, mat)
        
            % Efeito da duração da fixação no desempenho
            plotPSAfixDur(trlProps, mat)
        
            % Efeito da distância na performance (idealmente, invariante para pré-s mas variável para n-sac)
            % Talvez ordenar e usar quantis de modo que todos bins tenham mesma quantidade de pontos, e a distância média representada por eles?
            plotPSAprobesFixDist(trlProps, mat)
        end
    end

    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Gráficos de mapa de PSA');
        
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Gráficos de mapa de PSA');
            % Mapa de d-primes em função da distância e orientação, em relação a nSacc e probe
            plotPSAdprimeMap(probePos, probePosFix, nSaccProbePos, trlProps, mat)
        
            % Mapa de d-prime para verificar efeito da direção da sacada
            plotPSAdprimeMap(rotProbePos, rotProbeFix, rotNSaccProbePos, trlProps, mat, true);
        end
    end

    %% -- TAREFA DE FORRAGEAMENTO --
    fprintf('\n\n-- TAREFA DE FORRAGEAMENTO --\n')
    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Fixações na tarefa de forrageamento');
        
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Fixações na tarefa de forrageamento');
            % Distribuição das durações de fixação em ROI, comparada com a duração
            % das fixações propriamente
            plotFixDur(trlProps, mat, true)
        
            % Duração da fixação em função de ser alvo ou não e de ter acerto ou não
            plotForFixDur(trlProps, mat)

            % Frequência com que as sacadas são feitas para os estímulos mais
            % próximos ao pré-probe

        end
    end

    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Desempenho em tarefa de forrageamento');
        
        if strcmp(skipMode, 'all')
            skipAll = true;
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Desempenho em tarefa de forrageamento');
            % Efeito do número de fixações anteriores na duração da atual e no
            % desempenho do forrageamento
            plotForCountEffect(trlProps, mat)
            
            % Efeito de recência no desempenho do forrageamento
            % Ou quão isolado, faltaria fazer...
            plotForRecencyEffect(trlProps, mat)
        end
    end

    %% -- PROPRIEDADES DA TAREFA --
    disp('-- PROPRIEDADES DA TAREFA --')
    if ~skipAll
        [skipMode, closePlots] = promptPlotStep('Distância entre os estímulos');
        
        if strcmp(skipMode, 'all')
            skipAll = true; %#ok<NASGU> 
        elseif strcmp(skipMode, 'none')
            if closePlots, close all; end
            
            disp('> Distância entre os estímulos');
            %% Distância mínima entre estímulos
            taskStimDistr(trlProps, mat, 2)
        end
    end

%     if ~skipAll
%         [skipMode, closePlots] = promptPlotStep('??');
%         
%         if strcmp(skipMode, 'all')
%             skipAll = true;
%         elseif strcmp(skipMode, 'none')
%             if closePlots, close all; end
%             
%             disp('??');
%         end
%     end



end
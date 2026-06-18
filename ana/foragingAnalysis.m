function [trlProps, analysis, eyeData, evTimes] = foragingAnalysis(subj, ses)
    % Convenção:
    % vou chamar os estímulos forrageado (for), sacádico (sacc) e 
    % não-sacádico (nSacc) e identificá-los como -1, 0 e 1 com respeito ao 
    % tempo em que são vistos quanto ao início da sacada

    currFolder = fileparts(mfilename('fullpath')); parentFolder = fileparts(currFolder);
    addpath(genpath(fullfile(currFolder, 'dep'))); addpath(parentFolder);

    params = foragingParams;
    outFolder = fullfile(params.currFolder, params.outFolder);

    [mat, edf, sesStr] = foragingLoad(subj, ses, outFolder, params);

    mat.results.trialOrder(2,:,:) = mat.results.trialOrder(2,:,:) .* reshape(1:mat.tkP.nBlocks, 1, 1, []);

    %% Extrai as informações de cada trial, que já passaram por uma pré-seleção
    [trlProps, eyeData, evTimes] = foragingTrlProps(mat, edf, sesStr, subj);
    clear currFolder edf outFolder params parentFolder sesStr ses subj;

    trlProps = trlProps([trlProps.trlKeep]);

    fprintf('\n----------------------------------\nNúmero de trials pré-selecionados: %d\n', sum(logical([trlProps.trlKeep])))



    %% Rejeita trials com base no período de ruído rosa visto
    P3SaccLatencyLims = -1*[mat.prm.pinkNoiseDur+.2 .020];
    P3SaccLatency = [trlProps.saccInt1]/1000;
    goodTrl = P3SaccLatency >= P3SaccLatencyLims(1) & P3SaccLatency <= P3SaccLatencyLims(2);

    trlProps = trlProps(goodTrl);

    fprintf('Número final de trials selecionados: %d\n', sum(goodTrl))

    %% Análise exploratória do comportamento ocular
    fixPos        = pixel_to_dva([trlProps.preProbePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    fixPosFix     = pixel_to_dva([trlProps.preProbePosFixPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    probePos      = pixel_to_dva([trlProps.probePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    probePosFix   = pixel_to_dva([trlProps.probePosFixPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    nSaccProbePos = pixel_to_dva([trlProps.nSaccProbePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    foragingEyePlots(fixPos, fixPosFix, probePos, probePosFix, nSaccProbePos)




    %% Gráfico principal: Efeito pré-sacádico E tabela de contingência
    titlePSA  = 'Efeito pré-sacádico em tarefa de forrageamento';
    xlabelPSA = ["Forrageamento", "Sacádico", "Não-sacádico"];
    ylabelPSA = 'Acertos (%)';
    [main, counts] = getPSAeffect(trlProps);
    PSA.main = main; PSA.main.counts = counts;
    
    figure;
    b = bar(counts(1,:)./counts(2,:)*100, 'FaceColor', 'flat');
    b(1).CData(1,:) = [1 1 1];
    set(gca, 'TickDir', 'out', 'Box', 'off')
    xticklabels(xlabelPSA);
    ylabel(ylabelPSA); ylim([0 100]);
    title(titlePSA);


    %% Efeito da ordem das perguntas

    titleOrder  = 'Efeito de ordem das respostas';
    xlabelOrder = {["F", "S", "N"], ["S", "F", "N"], ["S", "N", "F"], ["S", "N"];
                   ["F", "N", "S"], ["N", "F", "S"], ["N", "S", "F"], ["N", "S"]};
    whiteIdx = [1 2 3; 1 2 3];
    ylabelOrder = {'Acertos (%)'; 'Acertos (%)'};
    [hits, counts] = getPSAOrder(trlProps);
    PSA.order.hits = hits; PSA.order.counts = counts;

    figure;
    for r = 1:2, for c = 1:4 %#ok<ALIGN> 
        subplotIdx = (r - 1) * 4 + c;
        subplot(2, 4, subplotIdx);
        
        b = bar(PSA.order.hits{r, c}, 'FaceColor', 'flat');
        
        set(gca, 'TickDir', 'out', 'Box', 'off');
        xticklabels(xlabelOrder{r, c});
        ylim([0 100]);
        
        if c == 1
            ylabel(ylabelOrder{r});
        end
        if c < 4
%             b.CData = repmat([0.4 0.4 0.4], 3, 1);
            b.CData(whiteIdx(r,c), :) = [1 1 1]; 
            b.EdgeColor = [0 0 0];
        end
    end; end
    sgtitle(titleOrder)


    %% Efeito de categorias
    getPSAcat([trlProps.preProbeCat], [trlProps.probeCat], [trlProps.probeHit], [trlProps.nSaccProbeHit]);

    %% Efeito do desempenho no forrageamento

    %% Efeito do número de vistos (no pré-s e na duração da fixação)

    %% Efeito da duração do ruído rosa (split)

    %% Efeito da duração da fixação no desempenho

    % E talvez gráficos mais descritivos:

    %% Amplitudes de sacada, 


end
function plotPSAforagingPerformance(trl, PSA, drP)
    correctFor = logical(PSA.main.for.idx{1,1} + PSA.main.for.idx{2,2});
    incorrectFor = ~correctFor;
    cFtrl = trl(correctFor);
    iFtrl = trl(incorrectFor);
    [~, cFcounts] = getPSAeffect(cFtrl);
    [~, iFcounts] = getPSAeffect(iFtrl);
    
    % d-prime robusto
    calcDPrime = @(hit, total) norminv(min(max(hit./total, 0.001), 0.999)) * sqrt(2);
    
    % Acurácia (%)
    plotDataAcc = cell(1, 2);
    plotDataAcc{1} = (cFcounts(1,2:3) ./ cFcounts(2,2:3)) * 100;
    plotDataAcc{2} = (iFcounts(1,2:3) ./ iFcounts(2,2:3)) * 100;
    
    % Sensibilidade (d')
    plotDataSens = cell(1, 2);
    plotDataSens{1} = calcDPrime(cFcounts(1,2:3), cFcounts(2,2:3));
    plotDataSens{2} = calcDPrime(iFcounts(1,2:3), iFcounts(2,2:3));
    
    % Contagem total de trialspara exibir no gráfico
    trialCounts = cell(1, 2);
    trialCounts{1} = cFcounts(2,2:3);
    trialCounts{2} = iFcounts(2,2:3);
    
    titleBase = 'Efeito pré-sacádico e desempenho no forrageamento';
    
    % 1. Gera figura de Acurácia
    renderPSAforagingPerformance(plotDataAcc, trialCounts, 'Acertos (%)', ...
        'PSA foraging effect - Accuracy', [titleBase, ' (Acurácia)'], drP, false);
        
    % 2. Gera figura de Sensibilidade
    renderPSAforagingPerformance(plotDataSens, trialCounts, 'Sensibilidade (d'')', ...
        'PSA foraging effect - Sensitivity', [titleBase, ' (Sensibilidade)'], drP, true);
end
function [pPSA, TPSA, permTPSA] = permAnalysisPSA(allTrlProps, allSubj, nPerm, params, doPlot)
    if nargin < 5, doPlot = 0; end
    % ----
    %    d-primes reais por sujeito -> 
    % -> diferenças de d-primes por sujeito ->
    % -> média da diferença
    % ----

    nSubj = numel(allSubj);
    dPerCond = zeros(nSubj, 2);

    %% Obtém a estatística para os dados
    for i = 1:nSubj
        subj = allSubj(i);
        subjTrlProps = allTrlProps([allTrlProps.subjNum] == subj);
        [main, ~] = getPSAeffect(subjTrlProps, 0);
        dPerCond(i,:) = [main.sacc.d, main.nSacc.d];
    end
    TPSA = sum(sum(dPerCond, 1).^2);

    %% Realiza permutações e obtém aa estatística para cada uma delas
    permTPSA = zeros(1,nPerm);
    for b = 1:nPerm
        permDPerCond = zeros(nSubj, 2);
        permAllTrlProps = foragingPSApermute(allTrlProps);
        for i = 1:nSubj
            subj = allSubj(i);
            subjTrlProps = permAllTrlProps([permAllTrlProps.subjNum] == subj);
            [main, ~] = getPSAeffect(subjTrlProps, 0);
            permDPerCond(i,:) = [main.sacc.d, main.nSacc.d];
        end
        permTPSA(b) = sum(sum(permDPerCond, 1).^2);
    end
    
    pPSA = plotPSAPermTest(TPSA, permTPSA, params, doPlot);
end
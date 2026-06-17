function [result, counts] = getPSAeffect(trlProps)

    [result.for.table, result.for.idx] = ...
        contTable([trlProps.forProbeCat], [trlProps.forProbeResp]);
    result.for.correct = result.for.table(1,1)+result.for.table(2,2);
    result.for.total   = sum(result.for.table(:));

    fprintf('\n\nTabela: condição pré-sacádica (pS)')
    [result.sacc.table, result.sacc.idx, result.sacc.d, result.sacc.c] = ...
        contTable([trlProps.probeCat], [trlProps.probeResp], 1);
    result.sacc.correct = trace(result.sacc.table);
    result.sacc.total   = sum(result.sacc.table(:));


    fprintf('\n\nTabela: condição não sacádica (nS)')
    [result.nSacc.table, result.nSacc.idx, result.nSacc.d, result.nSacc.c] = ...
        contTable([trlProps.nSaccProbeCat], [trlProps.nSaccProbeResp], 1);
    result.nSacc.correct = trace(result.nSacc.table);
    result.nSacc.total   = sum(result.nSacc.table(:));
    

    counts = [result.for.correct result.sacc.correct result.nSacc.correct; ...
        result.for.total result.sacc.total result.nSacc.total];
end
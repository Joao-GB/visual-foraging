function [result, counts] = getPSAeffect(trlProps)

    fprintf('\n\nTabela: forrageamento (F)')
    [result.for.table, result.for.idx, result.for.d, result.for.c] = ...
        contTable([trlProps.forProbeCat], [trlProps.forProbeResp], 1);
    result.for.correct = result.for.table(1,1)+result.for.table(2,2);
    result.for.total   = sum(result.for.table(:));

    fprintf('\n\nTabela: condição sacádica (S)')
    [result.sacc.table, result.sacc.idx, result.sacc.d, result.sacc.c] = ...
        contTable([trlProps.probeCat], [trlProps.probeResp], 1);
    result.sacc.correct = trace(result.sacc.table);
    result.sacc.total   = sum(result.sacc.table(:));


    fprintf('\n\nTabela: condição não sacádica (N)')
    [result.nSacc.table, result.nSacc.idx, result.nSacc.d, result.nSacc.c] = ...
        contTable([trlProps.nSaccProbeCat], [trlProps.nSaccProbeResp], 1);
    result.nSacc.correct = trace(result.nSacc.table);
    result.nSacc.total   = sum(result.nSacc.table(:));
    

    counts = [result.for.correct result.sacc.correct result.nSacc.correct; ...
        result.for.total result.sacc.total result.nSacc.total];
end
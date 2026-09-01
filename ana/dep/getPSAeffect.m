function [result, counts] = getPSAeffect(trlProps, showTable)
    if nargin < 2, showTable = true; end

    if showTable, fprintf('\n\nTabela: forrageamento (F)'); end
    [result.for.table, result.for.idx, result.for.d, result.for.c] = ...
        contTable([trlProps.forProbeCat], [trlProps.forProbeResp], showTable);
    result.for.correct = result.for.table(1,1)+result.for.table(2,2);
    result.for.total   = sum(result.for.table(:));

    if showTable, fprintf('\n\nTabela: condição sacádica (S)'); end
    [result.sacc.table, result.sacc.idx, result.sacc.d, result.sacc.c] = ...
        contTable([trlProps.probeCat], [trlProps.probeResp], showTable);
    result.sacc.correct = trace(result.sacc.table);
    result.sacc.total   = sum(result.sacc.table(:));


    if showTable, fprintf('\n\nTabela: condição não sacádica (N)'); end
    [result.nSacc.table, result.nSacc.idx, result.nSacc.d, result.nSacc.c] = ...
        contTable([trlProps.nSaccProbeCat], [trlProps.nSaccProbeResp], showTable);
    result.nSacc.correct = trace(result.nSacc.table);
    result.nSacc.total   = sum(result.nSacc.table(:));
    

    counts = [result.for.correct result.sacc.correct result.nSacc.correct; ...
        result.for.total result.sacc.total result.nSacc.total];
end
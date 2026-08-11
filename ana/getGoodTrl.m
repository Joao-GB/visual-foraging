function [allGoodTrl, keepIdx] =  getGoodTrl(trl, mat, trlName)
    % O próprio foragingTrlProps já rejeita os trials em que o probe sacádico 
    % é já visto por meio da linha:
    % trl(i).trlKeep = trlKeep && ~trl(i).probeSeen && hasProbeResp(i);
    keepIdx = logical([trl.trlKeep]);

    allGoodTrl = false(size(trl));
    
    % Calcula o limite de latência apenas para os trials mantidos
    P3SaccLatencyLims = -1 * [mat.prm.pinkNoiseDur + mat.prm.maxDelayFixOffP3, mat.prm.minP3Dur];
    P3SaccLatency = [trl(keepIdx).saccInt] / 1000;
    
    latencyMask = P3SaccLatency >= P3SaccLatencyLims(1) & P3SaccLatency <= P3SaccLatencyLims(2);
    
    % Atribui verdadeiro apenas onde AMBOS os critérios passaram
    allGoodTrl(keepIdx) = latencyMask;
    
    fprintf('\n----------------------------------\n%s\nNúmero de trials pré-selecionados: %d\n', trlName, sum(keepIdx))
    fprintf('Número de trials válidos final: %d\n', sum(allGoodTrl))
end
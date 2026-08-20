function [allGoodTrl, keepIdx] =  getGoodTrl(trl, mat, trlName)
    % O próprio foragingTrlProps já rejeita os trials em que o probe sacádico 
    % é já visto por meio da linha:
    % trl(i).trlKeep = trlKeep && ~trl(i).probeSeen && hasProbeResp(i);
    keepIdx = logical([trl.trlKeep]);

    allGoodTrl = false(size(trl));
    
    % Calcula o limite de latência apenas para os trials mantidos:
    % (i)  a diferença entre o início da sacada e o fim do estímulo não pode
    %     ser muito grande
    % (ii) o sujeito deve ver um mínimo do estímulo
    maxDelay = -mat.prm.maxDelayFixOffP3;
    minSeen  = - mat.prm.minP3Dur;
%     P3SaccLatencyLims = [mat.prm.maxDelayFixOffP3, mat.prm.minP3Dur];
    P3SaccInterval = [trl(keepIdx).saccInterval] / 1000;
    P3SaccLatency = [trl(keepIdx).saccLatency] / 1000;
    
    latencyMask = P3SaccInterval >= maxDelay & P3SaccLatency <= minSeen;
    
    % Atribui verdadeiro apenas onde AMBOS os critérios passaram
    allGoodTrl(keepIdx) = latencyMask;
    
    fprintf('\n----------------------------------\n%s\nNúmero de trials pré-selecionados: %d\n', trlName, sum(keepIdx))
    fprintf('Número de trials válidos final: %d\n', sum(allGoodTrl))
end
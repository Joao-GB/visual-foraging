function [allGoodTrl, keepIdx] =  getGoodTrl(trl, mat, trlName)
    % O próprio foragingTrlProps já rejeita os trials em que o probe sacádico 
    % é já visto por meio da linha:
    % trl(i).trlKeep = trlKeep && ~trl(i).probeSeen && hasProbeResp(i);
    keepIdx = logical([trl.trlKeep]);

    allGoodTrl = false(size(trl));
    
    %% 1. Calcula o limite de latência apenas para os trials mantidos:
    % (i)  a diferença entre o início da sacada e o fim do estímulo não pode
    %     ser muito grande
    % (ii) a diferença entre o fim da sacada (início da próxima fixação) e
    %     o fim da fase 3 tem que ser maior que um mínimo (uns 10 ou 15 ms)
    % (ii) o sujeito deve ver um mínimo do estímulo
    maxDelay = -mat.prm.maxDelayFixOffP3;
    minSeen  = - mat.prm.minP3Dur;
%     P3SaccLatencyLims = [mat.prm.maxDelayFixOffP3, mat.prm.minP3Dur];
    P3SaccInterval1= [trl(keepIdx).saccInterval1]/ 1000;
    P3SaccInterval = [trl(keepIdx).saccInterval] / 1000;
    P3SaccLatency = [trl(keepIdx).saccLatency] / 1000;
    
    % 
    latencyMask = P3SaccInterval >= maxDelay & P3SaccLatency <= minSeen & P3SaccInterval1 > .02; % & P3SaccInterval < 0

    %% 2. Identifica as distâncias entre fixação e probe
    maxDist = mat.prm.gaborSize_dva/2 + 1.5;

    probePos      = pixel_to_dva([trl(keepIdx).probePosPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';
    probePosFix   = pixel_to_dva([trl(keepIdx).probePosFixPix], 'dist', mat.prm.screenDist, 'width', mat.dpP.monitorW_mm/10, 'res', mat.dpP.screenRes.width)';

    distanceMask = vecnorm(probePos'-probePosFix') <= maxDist;
    
    % Atribui verdadeiro apenas onde AMBOS os critérios passaram
    allGoodTrl(keepIdx) = latencyMask & distanceMask;
    
    fprintf('\n----------------------------------\n%s\nNúmero de trials pré-selecionados: %d\n', trlName, sum(keepIdx))
    fprintf('Número de trials válidos final: %d\n', sum(allGoodTrl))
end
function allGoodTrl =  getGoodTrl(trl, mat)
    keepIdx = logical([trl.trlKeep]);

    allGoodTrl = false(size(trl));
    
    % Calcula o limite de latência apenas para os trials mantidos
    P3SaccLatencyLims = -1 * [mat.prm.pinkNoiseDur + 0.2, 0.020];
    P3SaccLatency = [trl(keepIdx).saccInt] / 1000; % Seguro: só puxa dos já filtrados
    
    latencyMask = P3SaccLatency >= P3SaccLatencyLims(1) & P3SaccLatency <= P3SaccLatencyLims(2);
    
    % Atribui verdadeiro apenas onde AMBOS os critérios passaram
    allGoodTrl(keepIdx) = latencyMask;
    
    % Mostra o print e filtra de uma vez só no final
    fprintf('\n----------------------------------\nNúmero de trials pré-selecionados: %d\n', sum(keepIdx))
    fprintf('Número de trials válidos final: %d\n', sum(allGoodTrl))
end
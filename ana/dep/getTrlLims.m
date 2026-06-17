function [trlLimIdx, trlLimAbs, trlDur] = getTrlLims(mat, messages, sesLimIdx, eventTimes)

    % Assumo que meu programa retorna fim para todos os trials, inclusive
    % se não for bom
    trlOnIdx = find(contains(messages(sesLimIdx(1):sesLimIdx(2)), mat.prm.msg.on.trl{2})) + sesLimIdx(1) - 1;
    trlOffIdx = zeros(1, length(trlOnIdx)); auxIdx = [trlOnIdx sesLimIdx(2)];
    for i=1:length(trlOnIdx)
        aux = find(contains(messages(auxIdx(i):auxIdx(i+1)), mat.prm.msg.off.trl{2}), 1) + trlOnIdx(i) - 1;
        if isempty(aux)  %% Apenas para qualquer erro que aconteça de eu não pausar direito 
            trlOffIdx(i) = trlOnIdx(i+1);
        else
            trlOffIdx(i) = aux;
        end
    end
    trlLimIdx = [trlOnIdx; trlOffIdx];
    trlLimAbs = reshape(eventTimes(1, trlLimIdx), size(trlLimIdx));
    trlDur = double(trlLimAbs(2,:)-trlLimAbs(1,:))/1000;
end
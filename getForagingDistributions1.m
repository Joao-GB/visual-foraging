function [nTs, nStims, targetOri, modTimes, nStimsToReport, orderToReportSets] = getForagingDistributions1(nStims, nTMin, nTMax, nTrials, nBlocks, params)
    % nStims = 18;

% (a) Distribuição da quantidade de alvos: binomial, de modo que (1) não requer
%     número de alvos; (2) em média, metade serão alvos; (3) um estímulo é 
%     alvo independentemente se o outro é; (4) não há como antecipar se o próximo a ser
%     olhado é alvo ou não acima da chance -- i.e., a predição ótima com base
%     no passado é sempre 1/2
        p = .5;
        nTs = binornd(nStims-1, p, nBlocks, nTrials);

    % (b) Distribuição da orientação-alvo: pseudo-uniforme nas orientações 
    %     de allOri, de modo que dois blocos consecutivos têm orientações
    %     distintas (claro, se houver mais de 2 alvos possíveis)
    aux = repmat(params.allOri, 1, ceil(nBlocks / length(params.allOri)));
    aux = aux(randperm(numel(aux)));
    if numel(params.allOri) > 2
        while(any(diff(aux) == 0)), aux = aux(randperm(numel(aux))); end
    end
    targetOri = aux(1:nBlocks);

    % (c) Distribuição do instante de modificação: após quantos alvos fixados 
    %     será apresentado o ruído rosa com orientação.
        modTimes = randsample(nTMin:nTMax, nBlocks*nTrials, true);
        modTimes = reshape(modTimes, nBlocks, nTrials);
%         modTimes = modTimes - 1;

    % (d) Distribuição da quantidade de estímulos cuja orientação deve ser
    %     reportada por trial: sempre pergunto de um não visto, às vezes de
    %     um já visto e os recém-fixados eu decido se pergunto na hora, com
    %     base na quantidade de estímulos fixados entre a fase 3 e PM.

        nStimsToReport = ones(params.maxToReport, nTrials, nBlocks);
        
        % Em cada trial, o sujeito deve reportar dos estímulos vistos
        % anteriormente, do atual ou dos não vistos em ordem aleatória
        [~, orderToReportSets] = sort(rand(3, nTrials, nBlocks), 1);
end
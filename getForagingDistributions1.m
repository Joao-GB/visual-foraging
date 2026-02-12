function [nTs, nStims, targetOri, modTimes, nStimsToReport, orderToReportSets] = getForagingDistributions1(nTMax, nTrials, nBlocks, params)
    nStims = 19;

% (a) Distribuição da quantidade de alvos: binomial, de modo que (1) não requer
%     número de alvos; (2) em média, metade serão alvos; (3) um estímulo é 
%     alvo independentemente se o outro é; (4) não há como antecipar se o próximo a ser
%     olhado é alvo ou não acima da chance -- i.e., a predição ótima com base
%     no passado é sempre 1/2
        p = .5;
        nTs = binornd(nStims, p, nBlocks, nTrials);

%         if nTargets > 1 && nTargets < nStims
%             numTargetsPMF = [.25 .50 .25];
%             targetsSpace  = nTargets + [-1 0 1];
%             nTs = randsample(targetsSpace, nBlocks*nTrials, true, numTargetsPMF);
%             nTs = reshape(nTs, [nBlocks nTrials]);
%         else
%             nTs = nTargets*ones(nBlocks, nTrials);
%         end

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
        
        % modeDistr = floor(3*nStims/4);
        % targetModTimePMF = robust_beta_pmf(nStims, modeDistr, 'peakness', 1.5);
        targetModTimePMF = modTime_pmf(nTMax);
        modTimes = randsample(1:nTMax, nBlocks*nTrials, true, targetModTimePMF);
        modTimes = reshape(modTimes, nBlocks, nTrials);
%         modTimes = modTimes - 1;

    % (d) Distribuição da quantidade de estímulos cuja orientação deve ser
    %     reportada por trial: sempre pergunto de um não visto, às vezes de
    %     um já visto e os recém-fixados eu decido se pergunto na hora, com
    %     base na quantidade de estímulos fixados entre a fase 3 e PM.

        nStimsToReport = round((params.maxToReport - params.minToReport)*rand(nBlocks,nTrials) + params.minToReport);
        nStimsPM = zeros(nBlocks, nTrials);
        nStimsPost = ones(nBlocks, nTrials);
        nStimsPre = ones(nBlocks, nTrials);

        % remainder = nStimsToReport - nStimsPM - nStimsPost;
        remainder = nStimsToReport - nStimsPM - nStimsPost - nStimsPre;

        % Garante que, se houver 2 ou mais estímulos no resto, pelo menos
        % um vá para o pós. Ao fim, espera-se que remainder seja de uns
        % nStimsPre = remainder >= 2;
        % remainder = remainder - nStimsPre;

        % Escolhe uniformemente quem do remainder vai para pré e pós
        addPre  = round(rand(size(remainder)).*remainder);
        addPost = remainder - addPre;

        nStimsPre = nStimsPre + addPre;
        nStimsPost = nStimsPost + addPost;

        nStimsToReport = cat(1, permute(nStimsPre,[3,2,1]), permute(nStimsPM, [3, 2, 1]), permute(nStimsPost,[3,2,1]));
        
        % Em cada trial, o sujeito deve reportar dos estímulos vistos
        % anteriormente, do atual ou dos não vistos em ordem aleatória
        [~, orderToReportSets] = sort(rand(3, nTrials, nBlocks), 1);
end
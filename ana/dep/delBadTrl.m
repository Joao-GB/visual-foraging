function [badTrl, nTrl, blkIdx] = delBadTrl(mat, messages, sesLimIdx, trlLimIdx, nTrl)

    % Encontra os trials descartados por pausa, timeout etc.
    badTrl = zeros(1, nTrl);
    for i=1:nTrl
        currTrlAllIdx = trlLimIdx(1,i):trlLimIdx(2,i);
        badCond = strcmp(messages(currTrlAllIdx), mat.prm.msg.err.trl{1}) | ... % Cobre efeitos de pausa
            ... strcmp(messages(currTrlIdx), mat.prm.msg.err.trl{2}) | ...   % Não incluo efeito de abort idx para manter nTrl 
            ...                                                              % por bloco. Essa info já está em trialOrder(2,:,:)
            ... strcmp(messages(currTrlIdx), mat.prm.msg.err.trl{3}) | ...   % Não é usado em foragingGabors o erro de NOP3
            strcmp(messages(currTrlAllIdx), mat.prm.msg.err.P1) | ...
            strcmp(messages(currTrlAllIdx), mat.prm.msg.err.P2);
        if any(badCond), badTrl(i) = true; end
    end

    % Encontra os blocos ruins, adelBadTrl saber
    % (i)  os que têm mensagem de erro dentro deles (sempre acompanhadas por 
    %     fim de bloco)
    % (ii) os que não têm mensagem de fim de bloco, e sempre estão
    %     associados a mensagens de interrupção de sessão
    % OBS: Veja que a análise não foi pensada para o caso de sessão recuperada
    blkOnIdx = find(contains(messages(sesLimIdx(1):sesLimIdx(2)), mat.prm.msg.on.blk{2})) + sesLimIdx(1) - 1;
    blkOffIdx = zeros(1, length(blkOnIdx)); nBlk = length(blkOnIdx);
    sesItMsg = sprintf(mat.prm.msg.off.ses{2}, mat.prm.msg.suffix{3});

    auxIdx = [blkOnIdx sesLimIdx(2)]; badBlk = zeros(1, nBlk);
    for i=1:nBlk
        currBlkIdx = auxIdx(i):auxIdx(i+1);
        blkOffIdx(i) = find(contains(messages(currBlkIdx), mat.prm.msg.off.blk{2})) + auxIdx(i) - 1;
        badCond = strcmp(messages(currBlkIdx), mat.prm.msg.err.blk) | ...    % (i)
                  strcmp(messages(currBlkIdx), sesItMsg);                    % (ii)
        if any(badCond), badBlk(i) = true; end
    end
    blkIdx = [blkOnIdx; blkOffIdx];
    badblkIdx = find(badBlk);
    for i=1:length(badblkIdx)
        badTrl(trlLimIdx(1,:) >= auxIdx(badblkIdx(i)) & trlLimIdx(1,:) <= auxIdx(badblkIdx(i)+1)) = true;
    end
    blkIdx(:, badblkIdx) = [];

    nTrl = sum(~badTrl);
    badTrl = logical(badTrl);
end
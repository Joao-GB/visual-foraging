function [pLI, pLA, pD] = getPhaseLims(m, st, tAI, mOn, mOff)
% Retorna os índices (pLI), tempo absoluto (pLA) e duração (pD) da fase
    pLI = [find(strcmp(m(tAI), mOn),1); find(strcmp(m(tAI), mOff),1)] + tAI(1) - 1;
    pLA = reshape(st(pLI), [2,1]);
    pD = double(pLA(2)-pLA(1))/1000;
end
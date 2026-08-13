function [T, tkP] = P3Onset5(tkP, prm, newFix)
% Versão para usar online, durante o trial
    tf = prm.pinkNoiseDur;
    medFix = median(tkP.fixQueue);
    % Se a mediana for muito curta, tirar esse valor faz o mínimo fixar
    % muito pequeno, por isso uso .8 como mínimo mínimo
    m = max(prm.minFixTime2, medFix - tf - prm.timeTosubtractFromMedianFix);
    M = medFix - tf;
    
    % Vou sortear o tempo uniformemente entre esse mínimo e a
    % mediana menos a duração do estímulo
    T = m + (M - m)*rand;
end
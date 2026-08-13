function [T, tkP] = P3Onset5OLD(tkP, prm, newFix)
% Versão para usar online, durante o trial
    tf = prm.pinkNoiseDur;
    medFix = median(tkP.fixQueue);
    auxMin = prm.minFixTimeBeforePink;
    m = min((medFix - tf), auxMin);
    M = max((medFix - tf), auxMin);
    
    % Em geral, eu espero que a mediana da fixação esteja bem distante do
    % mínimo que escolhi (algo em torno de 150 ms). Se for o caso, quer
    % dizer que vou sortear o tempo uniformemente entre esse mínimo e a
    % mediana menos a duração do estímulo; do contrário inverte quem é
    % mínimo e quem é máximo (então permitiria valores menores que o auxMin)
    T = m + (M - m)*rand;
end
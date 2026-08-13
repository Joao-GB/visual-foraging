function [T, tkP] = P3Onset4(tkP, prm, newFix)
% Versão para usar online, durante o trial
    tf = prm.pinkNoiseDur; 
    medFix = median(tkP.fixQueue);
    % estimativa do onset
    T = medFix - tf; % (D + tf);
end
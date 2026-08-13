function [T, tkP] = P3Onset2(tkP, prm, newFix)
% Versão para usar online, durante o trial
    beta  = prm.betaP3;

    tf = prm.pinkNoiseDur; 
    emaFix = tkP.fixProps.emaFix;

    emaFix = (1 - beta) * emaFix + beta * newFix;
    
    tkP.fixProps.emaFix = emaFix;

    % estimativa do onset
    T = emaFix - tf; % (D + tf);
end
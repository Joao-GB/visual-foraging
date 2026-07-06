function plotPSApinkDur(trl, mat)
    % Divisão dos trials com base no exactDur
    pinkNoiseDur = [trl.pinkNoiseDur];
    exactDur = floor(mat.prm.pinkNoiseDur * 1000) / 1000;

    figTitle = 'PSA pink noise duration';
    plotTitle = 'Efeito da duração do ruído rosa';

    plotPSAsplit(exactDur, pinkNoiseDur, trl, mat, figTitle, plotTitle)
end
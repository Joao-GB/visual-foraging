function pinkDurMask = plotPSApinkDur(trl, mat, figTitle, plotTitle)
    if nargin < 3, figTitle = 'PSA pink noise duration'; end
    if nargin < 4, plotTitle = 'Efeito da duração do ruído rosa'; end

    % Divisão dos trials com base no exactDur
    pinkNoiseDur = [trl.pinkNoiseDur];
    exactDur = floor(mat.prm.pinkNoiseDur * 1000) / 1000;

    pinkDurMask = plotPSAsplit(exactDur, pinkNoiseDur, trl, mat, figTitle, plotTitle);
end
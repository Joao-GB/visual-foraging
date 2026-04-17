function [fixCenter, stimCenters] = getStimLocations1(ROIparams, stimSize, minDist, ~, ~)
% O centro do hexágono não mais muda de posição, i.e., sempre há um
% estímulo no centro da tela
    totalRad = stimSize/2+2*minDist; %#ok<NASGU>

    hexCenter = [ROIparams(1); ROIparams(2)];%sampleCircleInRect(ROIparams, totalRad);
    hexAngle  = rand*60;

    stimCenters = getHexDistribution(hexCenter, minDist, 2, hexAngle)';

    nStims = size(stimCenters, 2);
    aux    = randsample(nStims,1);
    fixCenter = stimCenters(:,aux);
    stimCenters(:,aux) = [];
end
function [fixCenter, stimCenters] = getStimLocations1(ROIparams, stimSize, minDist, ~, ~)

    totalRad = stimSize/2+2*minDist;

    hexCenter = sampleCircleInRect(ROIparams, totalRad);
    hexAngle  = rand*60;

    stimCenters = getHexDistribution(hexCenter, minDist, 2, hexAngle)';

    nStims = size(stimCenters, 2);
    aux    = randsample(nStims,1);
    fixCenter = stimCenters(:,aux);
    stimCenters(:,aux) = [];
end
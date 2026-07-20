function aSigma = aSigmaFromStair(stair, prm)
    % aSigma = mean(stair.aSigma);
%     aux = min(stair.aSigma);
    aSigma = stair.aSigma - rem(stair.aSigma, prm.sigmaRem);
end
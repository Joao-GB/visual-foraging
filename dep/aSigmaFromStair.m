function aSigma = aSigmaFromStair(stair, prm)
    % aSigma = mean(stair.aSigma);
    aux = min(stair.aSigma);
    aSigma = aux - rem(aux, prm.sigmaRem);
end
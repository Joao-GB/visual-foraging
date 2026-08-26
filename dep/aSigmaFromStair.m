function aSigma = aSigmaFromStair(stair, prm, mode)
if nargin < 3
    mode = 0;
end
    % aSigma = mean(stair.aSigma);
%     aux = min(stair.aSigma);
if mode
    aSigma = stair - rem(stair, prm.sigmaRem);
else
    aSigma = stair.aSigma - rem(stair.aSigma, prm.sigmaRem);
end
end
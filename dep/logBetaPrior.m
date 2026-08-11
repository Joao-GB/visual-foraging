function [mu, sigma] = logBetaPrior(medianBeta, orders)
%LOGBETAPRIOR Parameters for a normal prior on log10(beta).
%
%   medianBeta : desired median of beta
%   orders     : orders of magnitude covered by +/- 1 SD
%
%   Example:
%       [mu, sigma] = logBetaPrior(0.5, 1);
%
%   gives a prior where the median beta is 0.5 and
%   +/- 1 SD corresponds to a factor of 10.

sigma = orders;
mu = log10(medianBeta);

end
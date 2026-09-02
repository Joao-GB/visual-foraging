function [Nmin, power] = binomial_power_N(acc, p0, alpha, targetPower, maxN)
% BINOMIAL_POWER_N
% Find the minimum number of trials required for a two-sided
% exact binomial test to achieve a desired statistical power.
%
% INPUTS
%   acc         = assumed true accuracy/proportion under H1
%   p0          = chance level (e.g. 0.5)
%   alpha       = significance level (e.g. 0.05)
%   targetPower = desired power (e.g. 0.80)
%   maxN        = maximum N to search (e.g. 10000)
%
% OUTPUTS
%   Nmin        = minimum number of trials achieving targetPower
%   power       = actual power at Nmin
%
% EXAMPLE
%   [N, power] = binomial_power_N(0.60, 0.50, 0.05, 0.80, 10000)

if nargin < 5
    maxN = 10000;
end

for N = 1:maxN

    % Probability of each possible number of correct responses
    k = 0:N;

    % Exact two-sided binomial test:
    % Determine which outcomes would be significant under H0.
    p_H0 = binopdf(k, N, p0);

    % Probability of each outcome under H0.
    % The exact binomial test rejects outcomes with sufficiently
    % small probability under H0.
    [p_sorted, idx] = sort(p_H0);

    cumulative = cumsum(p_sorted);
    significant = false(size(k));

    significant(idx(cumulative <= alpha)) = true;

    % Power = probability of significant outcomes under H1
    p_H1 = binopdf(k, N, acc);
    power = sum(p_H1(significant));

    if power >= targetPower
        Nmin = N;
        return
    end
end

error('Target power was not reached before maxN.');
end
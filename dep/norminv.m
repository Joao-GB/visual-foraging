function x = norminv(p, mu, sigma)

    if nargin < 2 || isempty(mu)
        mu = 0;
    end
    if nargin < 3 || isempty(sigma)
        sigma = 1;
    end

    x = mu + sigma .* sqrt(2) .* erfinv(2*p - 1);

end
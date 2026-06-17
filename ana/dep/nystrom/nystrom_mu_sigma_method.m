function [mu, sigma] = nystrom_mu_sigma_method(what, method)   
    if isempty(what)
        mu = NaN; sigma = NaN; return; 
    end
    if any(strcmp(method, {'mad', 'med', 'median'}))
        what(isnan(what)) = [];
        mu = median(what); sigma = 1.4826 * mad(what);
    elseif strcmp(method, 'mean')
        mu = mean(what, 'omitnan');   sigma = std(what, 'omitnan');
    end
end
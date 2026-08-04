function Nmin = minimumSubjects(K, delta, sigma, A, alpha, targetPower, doPlot)

    nVectors = sum([numel(K)>1, numel(delta)>1, numel(sigma)>1, numel(A)>1]);
    
    if nVectors > 1
        error('Exactly one of K, delta, sigma or A may be a vector.');
    end

    % Critical values
    zAlpha = norminv(1-alpha/2);
    zPower = norminv(targetPower);

    % Effective effect size
    deff = delta ./ sqrt(sigma.^2 + A./K);

    % Minimum number of subjects
    Nmin = ceil(((zAlpha + zPower)./deff).^2);

    % Plot if requested
    if doPlot

        % Determine which variable is varying
        if numel(K) > 1
            x = K;
            xlabelStr = 'Trials per subject (K)';
        elseif numel(delta) > 1
            x = delta;
            xlabelStr = '\delta';
        elseif numel(sigma) > 1
            x = sigma;
            xlabelStr = '\sigma';
        elseif numel(A) > 1
            x = A;
            xlabelStr = 'Measurement constant (A)';
        else
            error('Nothing to plot: all inputs are scalars.');
        end

        figure;
        plot(x,Nmin,'LineWidth',2)
        xlabel(xlabelStr)
        ylabel('Minimum number of subjects (N)')
        grid on

    end

end
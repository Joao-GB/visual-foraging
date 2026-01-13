function pmf = robust_beta_pmf(m, desired_mode, concentration_type, value)
    % concentration_type: 'variance', 'sum', or 'peakness'
    % value: corresponding parameter value
    
    p = (desired_mode - 1) / (m - 1);
    
    switch concentration_type
        case 'variance'
            % Direct variance control (0 to 0.25)
            target_var = min(max(value, 0.001), 0.2);
            sol = fsolve(@(x) beta_equations(x, p, target_var), [3, 3], optimset('Display','off'));
            alpha = max(sol(1), 1.01);
            beta = max(sol(2), 1.01);
            
        case 'sum'
            % Control α+β (higher = more peaked)
            k = max(value, 2.1);
            alpha = p * (k - 2) + 1;
            beta = k - alpha;
            
        case 'peakness'
            % Simple intuitive control (1-10 scale)
            k = 2 + value * 8;  % Map to [2, 10]
            alpha = p * (k - 2) + 1;
            beta = k - alpha;
    end
    
    % Generate distribution
    x_continuous = (0.5:(m-0.5)) / m;
    densities = betapdf(x_continuous, alpha, beta);
    pmf = densities / sum(densities);
end
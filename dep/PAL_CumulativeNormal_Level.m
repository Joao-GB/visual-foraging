function y = PAL_CumulativeNormal_Level(params, x, targetLevel, varargin)
    % PAL_CumulativeNormal_Level  Cumulative Normal PF parameterized by a target level.
    %
    % params: [alpha, beta, gamma, lambda]
    %         alpha now represents the stimulus intensity at 'targetLevel'
    % x: stimulus intensity (or probability, if 'Inverse' is called)
    % targetLevel: The desired overall probability (e.g., 0.85)
    
    [alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);
    
    % 1. Calculate the shifted z-score needed to anchor alpha at the targetLevel
    % c_target is the sensory proportion required to hit the targetLevel
    c_target = (targetLevel - gamma) ./ (1 - gamma - lambda);
    
    % Ensure target is mathematically reachable given guess and lapse rates
    if any(c_target(:) <= 0 | c_target(:) >= 1)
        error('targetLevel must be strictly between gamma and 1 - lambda');
    end
    
    % z_target is the standard normal z-score corresponding to c_target
    z_target = -sqrt(2) .* erfinv(1 - 2 .* c_target);
    
    % Z is the shifted input. When x == alpha, Z == z_target.
    Z = beta .* (x - alpha) + z_target;
    
    if ~isempty(varargin)
        if strncmpi(varargin{1}, 'Inverse', 3)
            % 'x' is actually the target probability 'y' here
            c = (x - gamma) ./ (1 - gamma - lambda);
            
            % Solve for stimulus intensity: Z = -sqrt(2) * erfinv(1 - 2c)
            y = alpha - (sqrt(2) .* erfinv(1 - 2 .* c) + z_target) ./ beta;
        end
        if strncmpi(varargin{1}, 'Derivative', 3)
            % Derivative of the PF with respect to x (used by Psi to find optimal placement)
            % d/dx [Phi(beta*(x-alpha) + z_target)] = beta * phi(...)
            y = (1 - gamma - lambda) .* beta .* (1/sqrt(2*pi)) .* exp(-(Z.^2) ./ 2);
        end
    else
        % Forward PF evaluation
        y = gamma + (1 - gamma - lambda) .* .5 .* erfc(-Z ./ sqrt(2));
    end
end
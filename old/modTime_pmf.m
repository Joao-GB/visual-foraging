function p = modTime_pmf(N, alpha)
% alpha indica quanto por cento da distribuição uniforme fica nas pré-pontas
if nargin < 2, alpha = .5; end
    if N < 3, p = ones(1,N)/N; return; end

    p = zeros(1,N);
    if N <= 4
        p(2:N-1) = 1/(N-2);
        return
    end
    M = N-4;
    p_int = 1/(M+2*alpha);
    p(3:N-2) = p_int;
    p(2) = alpha*p_int;
    p(N-1) = alpha*p_int;
end
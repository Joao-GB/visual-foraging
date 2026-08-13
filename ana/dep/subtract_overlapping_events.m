% Remove de A os eventos que se sobrepõem com algum evento de B.
% INPUTS:    
% A:         conjunto do qual se removerá eventos.
% B:         conjunto dos eventos a serem procurados em A.
% mode:      0 se não permitir sobreposição de fronteiras, 1 se permitir
%            (padrão: 0).
% OUTPUTS: 
%  C:        subconjunto de A dos eventos que não se sobrepõem com os de B.
% A_minus_C: subconjunto de A dos eventos que se sobrepõem com os de B, de
%            modo que C \cup A_minus_C = A (a menos da ordem dos eventos).

% Exemplos:
% A = [1, 5, 12; 4, 8, 15];
% B = [2, 7; 3, 12];
% subtract_overlapping_events(A, B, 0) --> retorna []
% subtract_overlapping_events(A, B, 1) --> retorna [12; 15]

function [C, A_minus_C] = subtract_overlapping_events(A, B, mode)
    if nargin < 3
        mode = 0;
    end
    % Initialize output
    C = [];
    A_minus_C = [];
    
    if isempty(A) || isempty(B)
        C = A;
        A_minus_C = [];
        return
    end
    
    % For each A event, check if it overlaps any B event
    for a_idx = 1:size(A, 2)
        a_start = A(1, a_idx);
        a_end = A(2, a_idx);
        
        % Check if A overlaps any B
        if mode == 0
            overlaps = any((B(1, :) <= a_end) & (B(2, :) >= a_start));
        else
            overlaps = any((B(1, :) < a_end) & (B(2, :) > a_start));
        end
        % If no overlap, keep the A event
        if ~overlaps
            C = [C, [a_start; a_end]];
        else
            A_minus_C = [A_minus_C, [a_start; a_end]];
        end
    end
end
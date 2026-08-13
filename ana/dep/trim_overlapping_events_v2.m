% Difere do trim_overlapping_events pois permite arrays com 3 linhas, para
% índices
% Faz com que os eventos do conjunto A terminem antes ou comecem depois de
% modo que não se sobreponham temporalmente com os eventos do conjunto B,
% retornando A modificado. Eventos definidos em tempo discreto (naturais).
% 
% Exemplo de input e output para entender como funciona: 1a linha é início
% dos eventos e 2a linha o fim.
% INPUTS:
% A = [1, 5, 10;
%      4, 8, 15];
% B = [2, 7;
%      3, 11];
% OUTPUT:
% C = [1, 4, 5, 12;
%      1, 4, 6, 15]
%
% Disclaimer: escrito pelo Deepseek, mas pensado em detalhe, quanto aos 
%             inputs e outputs desejados, por mim
function C = trim_overlapping_events_v2(A, B, mode)
    % Set default mode to 0 if not provided
    if nargin < 3
        mode = 0;
    end
    
    if isempty(B)
        C = A;
        return;
    end
    
    % Sort intervals
    [~, idx] = sort(A(1, :));
    A = A(:, idx);
    [~, idx] = sort(B(1, :));
    B = B(:, idx);
        
    rows = size(A, 1);
    % Merge overlapping B intervals
    B_merged = merge_intervals(B);
    C = [];
    
    for a_idx = 1:size(A, 2)
        a_start = A(1, a_idx);
        a_end = A(2, a_idx);
        if rows == 3
            a_type = A(3, a_idx);
        else
            a_type = NaN;
        end
            
        % Find B intervals that overlap with current A
        overlapping_Bs = B_merged(:, B_merged(1, :) <= a_end & B_merged(2, :) >= a_start);
        
        if isempty(overlapping_Bs)
            C = [C, aux_fun(a_start, a_end, a_type)];
            continue;
        end
        
        current_start = a_start;
        for b_idx = 1:size(overlapping_Bs, 2)
            b_start = overlapping_Bs(1, b_idx);
            b_end = overlapping_Bs(2, b_idx);
            
            % Mode 0: Trim overlaps completely (original behavior)
            if mode == 0
                if current_start < b_start
                    C = [C, aux_fun(current_start, b_start - 1, a_type)];
                end
                current_start = b_end + 1;
            
            % Mode 1: Keep boundaries (new behavior)
            else
                if current_start < b_start
                    C = [C, aux_fun(current_start, b_start, a_type)];
                end
                current_start = b_end;
            end
        end
        
        % Add the remaining segment after last B (if any)
        if current_start <= a_end
                C = [C, aux_fun(current_start, a_end, a_type)];
        end
    end
end
    
function vec = aux_fun(a1, a2, a3)
    if isnan(a3)
        vec = [a1; a2];
    else
        vec = [a1; a2; a3];
    end
end

function merged = merge_intervals(intervals)
    if isempty(intervals)
        merged = [];
        return;
    end
    
    [~, idx] = sort(intervals(1, :));
    intervals = intervals(:, idx);
    merged = intervals(:, 1);
    
    for i = 2:size(intervals, 2)
        last_end = merged(2, end);
        current_start = intervals(1, i);
        current_end = intervals(2, i);
        
        if current_start <= last_end + 1
            merged(2, end) = max(last_end, current_end);
        else
            merged = [merged, [current_start; current_end]];
        end
    end
end
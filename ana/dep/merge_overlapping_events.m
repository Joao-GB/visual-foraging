% Une os eventos de A com alguma intersecção
function B = merge_overlapping_events(A)
    if isempty(A)
        B = A;
        return
    end
    % Step 1: Sort intervals by start time
    sorted_A = reorder_events(A);
    
    % Step 2: Initialize merged list
    merged = [];
    current_start = sorted_A(1, 1);
    current_end = sorted_A(2, 1);
    
    % Step 3: Iterate and merge
    for i = 2:size(sorted_A, 2)
        next_start = sorted_A(1, i);
        next_end = sorted_A(2, i);
        
        % Check overlap or adjacency
        if next_start <= current_end
            % Merge: extend current_end if needed
            current_end = max(current_end, next_end);
        else
            % No overlap: save current interval and reset
            merged = [merged, [current_start; current_end]];
            current_start = next_start;
            current_end = next_end;
        end
    end
    
    % Add the last merged interval
    B = [merged, [current_start; current_end]];
end
function A = reorder_events(A)
    if isempty(A); return; end
    [starts_sorted, sort_idx] = sort(A(1, :));
    L = size(A,1);
    A(1, :) = starts_sorted;
    if L > 1
        for i=2:L
            row_sorted = A(i, sort_idx);
            A(i,:) = row_sorted;
        end
    end
end
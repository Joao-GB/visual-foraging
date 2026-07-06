
function [H, C] = getPSAorderOLD(trl)
    orders = reshape([trl.allProbesOrder]', 3, [])';
    hits   = reshape([trl.allHit]', 3, [])';

    % Reference permutations mapping to your target grid layout
    % Row 1 targets: [-1 0 1], [0 -1  1], [0 1 -1]
    % Row 2 targets: [-1 1 0], [1 -1  0], [1 0 -1]
    allPerms = [ ...
        -1  0  1;   % Perm 1
         0 -1  1;   % Perm 2
         0  1 -1;   % Perm 3
        -1  1  0;   % Perm 4
         1 -1  0;   % Perm 5
         1  0 -1;   % Perm 6
    ];

    % Initialize output cells as 2x4 grids
    H = cell(2,4);
    C = cell(2,4);
    
    % Base-3 key trick to group the 6 main permutations
    keyWeight = [9, 3, 1]; 
    rowKeys   = (orders + 2) * keyWeight';
    permKeys  = (allPerms + 2) * keyWeight';
    
    % --- POPULATE COLUMNS 1, 2, AND 3 ---
    % Loop through our 6 explicitly ordered permutations
    for i = 1:6
        % Determine grid coordinate based on permutation index
        if i <= 3
            r = 1; c_idx = i;     % Permutations 1, 2, 3 -> Grid Row 1
        else
            r = 2; c_idx = i - 3; % Permutations 4, 5, 6 -> Grid Row 2
        end
        
        idx = find(rowKeys == permKeys(i));
        N = numel(idx);
        
        if N > 0
            % sum(..., 1) ensures it stays a 1x3 vector even if only 1 trial matches
            C{r, c_idx} = sum(hits(idx, :), 1);
            H{r, c_idx} = (C{r, c_idx} / N) * 100; % Percentage
        else
            C{r, c_idx} = [0, 0, 0];
            H{r, c_idx} = [0, 0, 0];
        end
    end
    
    % --- POPULATE COLUMN 4 (Category Relationships) ---
    % Grid Row 1, Col 4: Category 0 comes before 1 (ignoring -1)
    % This matches Permutations 1, 2, and 3 completely!
    idx_0_before_1 = find(rowKeys == permKeys(1) | rowKeys == permKeys(2) | rowKeys == permKeys(3));
    [H{1,4}, C{1,4}] = getCatPairHits(orders(idx_0_before_1,:), hits(idx_0_before_1,:), 0, 1);

    % Grid Row 2, Col 4: Category 1 comes before 0 (ignoring -1)
    % This matches Permutations 4, 5, and 6 completely!
    idx_1_before_0 = find(rowKeys == permKeys(4) | rowKeys == permKeys(5) | rowKeys == permKeys(6));
    [H{2,4}, C{2,4}] = getCatPairHits(orders(idx_1_before_0,:), hits(idx_1_before_0,:), 1, 0);
end

% --- Helper Function to compute hits for a specific category pair sequence ---
function [pct, cnt] = getCatPairHits(subOrders, subHits, firstCat, secondCat)
    if isempty(subOrders)
        cnt = [0, 0]; pct = [0, 0]; return;
    end
    
    % Find column positions (1, 2, or 3) for both categories across all rows
    [~, posFirst]  = max(subOrders == firstCat, [], 2);
    [~, posSecond] = max(subOrders == secondCat, [], 2);
    
    % Linear indexing to extract the exact hit state for those positions
    rowIdx = (1:size(subHits, 1))';
    hitsFirst  = subHits(sub2ind(size(subHits), rowIdx, posFirst));
    hitsSecond = subHits(sub2ind(size(subHits), rowIdx, posSecond));
    
    % Total hits and percentage array: [FirstCategory, SecondCategory]
    cnt = [sum(hitsFirst), sum(hitsSecond)];
    pct = (cnt / size(subHits, 1))*100;
end
function [H, C, D] = getPSAorder(trl)
    orders = reshape([trl.allProbesOrder]', 3, [])';

    % Reference permutations mapping to your target grid layout
    % Row 1 targets: [-1 0 1], [0 -1 1], [0 1 -1]
    % Row 2 targets: [-1 1 0], [1 -1 0], [1 0 -1]
    allPerms = [ ...
        -1  0  1;   % Perm 1
         0 -1  1;   % Perm 2
         0  1 -1;   % Perm 3
        -1  1  0;   % Perm 4
         1 -1  0;   % Perm 5
         1  0 -1;   % Perm 6
    ];
    permutations = [ ...
         1 2 3;
         2 1 3;
         2 3 1;
         1 3 2;
         3 1 2;
         3 2 1;
    ];

    permsAux = [ ...
        1 2;
        2 1;
    ];

    % Initialize output cells as 2x4 grids
    H = cell(2,4);
    C = cell(2,4);
    D = cell(2,4);

    for i=1:6

        if i <= 3
            r = 1; c = i;     % Permutations 1, 2, 3 -> Grid Row 1
        else
            r = 2; c = i - 3; % Permutations 4, 5, 6 -> Grid Row 2
        end
        currOrderIdx = ismember(orders, allPerms(i,:), 'rows');
        if isempty(currOrderIdx)
            C{r,c} = [0 0 0];
            H{r,c} = [0 0 0];
            D{r,c} = [0 0 0];
        else
            [currResult, currCounts] = getPSAeffect(trl(currOrderIdx));
            aux = currCounts(1,:);
            C{r,c} = aux(permutations(i,:));
            aux = [currResult.for.correct currResult.sacc.correct currResult.nSacc.correct] ./ [currResult.for.total currResult.sacc.total currResult.nSacc.total]*100;
            H{r,c} = aux(permutations(i,:));
            aux = [currResult.for.d currResult.sacc.d currResult.nSacc.d];
            D{r,c} = aux(permutations(i,:));
        end

    end

    SNorderIdx = ismember(orders, allPerms(1,:), 'rows') | ismember(orders, allPerms(2,:), 'rows') | ismember(orders, allPerms(3,:), 'rows');
    if isempty(SNorderIdx)
        C{1,4} = [0 0]; H{1,4} = [0 0]; D{1,4} = [0 0];
    else
        [currResult, currCounts] = getPSAeffect(trl(SNorderIdx));
        aux = currCounts(1,2:3);
        C{1,4} = aux(permsAux(1,:));
        aux = [currResult.sacc.correct currResult.nSacc.correct] ./ [currResult.sacc.total currResult.nSacc.total]*100; 
        H{1,4} = aux(permsAux(1,:));
        aux = [currResult.sacc.d currResult.nSacc.d];
        D{1,4} = aux(permsAux(1,:));
    end

    NSorderIdx = ismember(orders, allPerms(4,:), 'rows') | ismember(orders, allPerms(5,:), 'rows') | ismember(orders, allPerms(6,:), 'rows');
    if isempty(NSorderIdx)
        C{2,4} = [0 0]; H{2,4} = [0 0]; D{2,4} = [0 0];
    else
        [currResult, currCounts] = getPSAeffect(trl(NSorderIdx));
        aux = currCounts(1,2:3);
        C{2,4} = aux(permsAux(2,:));
        aux = [currResult.sacc.correct currResult.nSacc.correct] ./ [currResult.sacc.total currResult.nSacc.total]*100; 
        H{2,4} = aux(permsAux(2,:));
        aux = [currResult.sacc.d currResult.nSacc.d];
        D{2,4} = aux(permsAux(2,:));
    end
end
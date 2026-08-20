function outLims = biggest_intersection_events(aLims, bLims)

    % Compute overlap with each interval in bLims
    starts = max(aLims(1), bLims(1,:));
    ends   = min(aLims(2), bLims(2,:));

    overlap = max(0, ends - starts);

    [m, idx] = max(overlap);

    if m == 0
        outLims = [];      % no positive overlap
    else
        outLims = [starts(idx); ends(idx)];
    end
end
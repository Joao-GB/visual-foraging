function pairs = pairStartEnd(starts, ends, eventLims)
    eventStart = eventLims(1);
    eventEnd   = eventLims(2);
    
    % Keep only indices within the event
    starts = starts(starts >= eventStart & starts <= eventEnd);
    ends   = ends(ends >= eventStart & ends <= eventEnd);
    
    nStarts = numel(starts);
    
    pairs = nan(nStarts,2);
    pairs(:,1) = starts(:);
    
    endIdx = 1;
    
    for k = 1:nStarts
    
        % Upper bound for a valid end
        if k < nStarts
            nextStart = starts(k+1);
        else
            nextStart = inf;
        end
    
        % Check whether the next unused end belongs to this start
        if endIdx <= numel(ends) && ...
           ends(endIdx) > starts(k) && ...
           ends(endIdx) < nextStart
    
            pairs(k,2) = ends(endIdx);
            endIdx = endIdx + 1;
    
        else
            % No matching end
            if k < nStarts
                pairs(k,2) = nextStart - 1;
            else
                pairs(k,2) = eventEnd;
            end
        end
    
    end
    pairs = pairs';
end
function split_array = split_events(indices, split_points, mode)
% SPLIT_EVENTS Split an array at specified points
%   split_array = split_events(indices, split_points) splits the indices
%   array at the specified split_points, placing split points at the
%   beginning of each new subarray (default mode)
%
%   split_array = split_events(indices, split_points, mode) allows
%   specifying the placement of split points:
%       'start' - split points at start of subarrays (default)
%       'end'   - split points at end of subarrays
%
%   Example:
%       indices = [4 5 6 7 8 9 10 11];
%       split_points = [7, 11];
%       
%       % Default mode ('start')
%       result = split_events(indices, split_points)
%       % Returns: {[4 5 6], [7 8 9 10], [11]}
%       
%       % 'end' mode
%       result = split_events(indices, split_points, 'end')
%       % Returns: {[4 5 6 7], [8 9 10 11]}
%
% Disclaimer: escrito pelo Deepseek, mas pensado em detalhe, quanto aos 
%             inputs e outputs desejados, por mim (diretor executivo)

    % Set default mode if not provided
    if nargin < 3
        mode = 'start';
    end
    
    % Validate inputs
    if ~isnumeric(indices) || ~isvector(indices)
        error('indices must be a numeric vector');
    end
    
    if ~isnumeric(split_points) || ~isvector(split_points)
        error('split_points must be a numeric vector');
    end
    
    % Convert to row vectors for consistency
    indices = indices(:)';
    split_points = split_points(:)';
    
    % Remove any split points not in the indices array
    split_points = intersect(split_points, indices);
    
    if isempty(split_points)
        split_array = {indices};
        return;
    end
    
    switch lower(mode)
        case 'start'
            split_array = split_at_start(indices, split_points);
        case 'end'
            split_array = split_at_end(indices, split_points);
        otherwise
            error('Mode must be either ''start'' or ''end''');
    end
end

function split_array = split_at_start(indices, split_points)
    % Split with split points at the start of each subarray
    
    % Find positions of split points
    [~, split_positions] = ismember(split_points, indices);
    split_positions = sort(split_positions); % Ensure they're in order
    
    % Add start and end boundaries
    boundaries = [1, split_positions, length(indices) + 1];
    
    split_array = cell(1, length(boundaries) - 1);
    
    for i = 1:length(boundaries) - 1
        start_idx = boundaries(i);
        end_idx = boundaries(i + 1) - 1;
        
        % Only include non-empty segments
        if start_idx <= end_idx
            split_array{i} = indices(start_idx:end_idx);
        end
    end
    
    % Remove any empty cells (shouldn't happen with proper inputs)
    split_array = split_array(~cellfun('isempty', split_array));
end

function split_array = split_at_end(indices, split_points)
    % Split with split points at the end of each subarray
    
    % Find positions of split points
    [~, split_positions] = ismember(split_points, indices);
    split_positions = sort(split_positions); % Ensure they're in order
    
    % Add start and end boundaries
    boundaries = [0, split_positions, length(indices)];
    
    split_array = cell(1, length(boundaries) - 1);
    
    for i = 1:length(boundaries) - 1
        start_idx = boundaries(i) + 1;
        end_idx = boundaries(i + 1);
        
        % Only include non-empty segments
        if start_idx <= end_idx
            split_array{i} = indices(start_idx:end_idx);
        end
    end
    
    % Remove any empty cells (shouldn't happen with proper inputs)
    split_array = split_array(~cellfun('isempty', split_array));
end
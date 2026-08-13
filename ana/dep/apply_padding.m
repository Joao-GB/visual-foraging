function X_out = apply_padding(X, pad_len, mode, padding_type, pad_value)
% APPLY_PADDING Add or remove padding to/from each row of a matrix.
% 
%   X_out = apply_padding(X, pad_len, mode, padding_type, pad_value)
%
% Inputs:
%   X            - Input matrix (N x T), where N is the number of signals
%   pad_len      - Number of points to pad/remove on each side
%   mode         - 'add' or 'remove'
%   padding_type - 'replicate', 'symmetric', or 'constant'
%   pad_value    - (Optional) scalar, only used if padding_type is 'constant'
%
% Output:
%   X_out        - Output matrix with padding added or removed

if nargin < 5
    pad_value = 0; % Default for constant padding
end

[N, T] = size(X);

switch lower(mode)
    case 'add'
        X_out = zeros(N, T + 2*pad_len);
        for i = 1:N
            row = X(i,:);
            switch lower(padding_type)
                case 'replicate'
                    left = repmat(row(1), 1, pad_len);
                    right = repmat(row(end), 1, pad_len);
                case 'symmetric'
                    left = row(pad_len+1:-1:2);
                    right = row(end-1:-1:end-pad_len);
                case 'constant'
                    left = pad_value * ones(1, pad_len);
                    right = pad_value * ones(1, pad_len);
                otherwise
                    error('Unknown padding_type');
            end
            X_out(i,:) = [left, row, right];
        end
        
    case 'remove'
        if 2*pad_len >= T
            error('Padding length too large to remove from data.');
        end
        X_out = X(:, pad_len+1:end-pad_len);

    otherwise
        error('Mode must be ''add'' or ''remove''');
end

end

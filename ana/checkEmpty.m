function out = checkEmpty(a, m)
if nargin < 2, m = 0; end
 if m
     out = a.*2-1;
 else
    if isempty(a), out = 0; 
    else, out = a;
    end
end
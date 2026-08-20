function out = diff_events(A)

if    isempty(A), out = 0;
else,             out = diff(A); end
end
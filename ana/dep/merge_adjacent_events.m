% Faz com que eventos do conjunto A que são avizinhados (segundo limiar thr)
% virem um só.
% ASSUME QUE NÃO há sobreposição, por isso deve fazer mais sentido usar
% após merge_overlapping_events
function [A] = merge_adjacent_events(A_on, A_off, thr)
if isempty(A_on)
    A = []; return
end
if nargin < 3
    thr   = A_off;
    A_off = A_on(2,:);
    A_on  = A_on(1,:);
end
cnt = 1;
while cnt > 0
%     disp('while...')
    b = mov_props(A_on, A_off);
    itv_bl = length(b.interval);
    
    cnt = 0;
    for i = 1:itv_bl
        if b.interval(i) <= thr
%             A_on(i) = min(A_on(i), A_on(i+1));
%             A_off(i+1) = max(A_off(i), A_off(i+1));
            A_off(i)  = NaN;
            A_on(i+1) = NaN;
            cnt = cnt + 1;
        end
    end
    A_off(isnan(A_off)) = [];
    A_on(isnan(A_on))   = [];
    A = [A_on; A_off];
end
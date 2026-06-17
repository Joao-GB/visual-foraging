function [fixDur, fixLims] = foragingGetFix(stmOn, stmOff, ev, evT)
    % Essa versão usa informações de início e fim de sacadas do próprio
    % Eyelink
    L = length(stmOn);
    fixDur = zeros(1, L); fixLims = cell(1, L);
    for i=1:L
        currEv = ev(stmOn(i):stmOff(i));
        stFix = find(cellfun(@(x) strcmp('STARTFIX', x), currEv)) + stmOn(i) - 1;
        enFix = find(cellfun(@(x) strcmp('ENDFIX', x), currEv)) + stmOn(i) - 1;
        if isempty(stFix)
            aux = find(cellfun(@(x) strcmp('STARTFIX', x), ev(1:stmOn(i))), 1, 'last');
            stFix = [stFix aux]; %#ok<AGROW>
        end
        if isempty(enFix) || stFix(end) > enFix(end)
            aux = find(cellfun(@(x) strcmp('ENDFIX', x), ev(stmOff(i):end)), 1, 'first') + stmOff(i) - 1;
            enFix = [enFix aux]; %#ok<AGROW>
        end
%         if length(enFix) > 1, disp('Microssacadas?'); end

        fixLims{i} = [evT(1,enFix);evT(2,enFix)];
        fixDur(i) = double(evT(2,enFix(end)) - evT(1,stFix(1)))/1000;
    end
end
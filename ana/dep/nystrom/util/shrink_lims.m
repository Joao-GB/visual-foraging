function [new_lims, rm] = shrink_lims(lims)
    if isempty(lims)
        new_lims = [];
        rm = [];
        return;
    end
    l_on = lims(1,:);
    l_off = lims(2,:);
    l_on = l_on + 1;
    l_off = l_off - 1;
    aux = l_off - l_on;
    rm = lims(:, aux(aux < 0));
    new_lims = [l_on; l_off];
end
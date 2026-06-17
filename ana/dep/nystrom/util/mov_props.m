% Obtém as durações, intervalos e índices associados a um evento em tempo
% discreto, dados seus onset e offset
function m = mov_props(on, off)
    % Posso dar uma entrada com formato [on; off]. Isso dá flexibilidade ao
    % programa enquanto eu não altero as demais funções que a utilizam, mas
    % com on e off separados...
    if isempty(on)
        m.lims     = [];
        m.duration = [];
        m.interval = [];
        m.idx      = [];
        return
    end
    if ~exist('off', 'var')
        off = on(2,:);
        on  = on(1,:);
    end
    % Remove NaNs
    on(isnan(on)) = []; off(isnan(off)) = [];
    
    m.lims = [on; off];
    m.duration = off - on + 1;
    L = length(on);
    if L >= 2
        m.interval = on(2:end) - off(1:end-1);
    else
        m.interval = NaN;
    end

    m.idx = [];
    for i = 1:L
        m.idx = [m.idx m.lims(1, i):m.lims(2, i)];
    end

    if exist('d', 'var')
        L = length(d);
        m.data = d;
        m.data(:, setdiff(1:L, m.idx)) = NaN;
    end
end
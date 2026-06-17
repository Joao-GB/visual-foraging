% Obtém as piscadas/ruídos.
% Talvez aprimorar detecção com dados de pupila, apesar de os ruídos
% aparentemente se sobreporem...
function b = get_blinks(data, par, thr, f)
% INPUT:
% data      = deve estar em dva, para que o ruído esteja limitado por 90

    if ~exist('f','var')
        f = 20;
    end

    % Se limiar for radial
    if length(thr) == 1
        thr = [thr thr];
    end
    
    x = double(data(1, :));
    y = double(data(2, :));

    L = length(x);

    % 1) Encontro os índices em que ao menos uma das coordenadas excede o limiar
    aux(1, :) = abs(x) > thr(1);
    aux(2, :) = abs(y) > thr(2);
    aux = aux(1, :) | aux(2, :);

    [b_on, b_off] = get_start_stop(aux);
    x_filt = medfilt1(x, f);
    y_filt = medfilt1(y, f);

    % 2) Encontra os mínimos locais avizinhados do início e do fim da
    % piscada, mas nos dados filtrados para garantir que mínimos locais
    % representem atividade biológica
    for i=1:length(b_on)
        aux_x = b_on(i);
        aux_y = b_on(i);
        while aux_x >= 0 && x_filt(aux_x) - x_filt(aux_x - 1) >= 0
            aux_x = aux_x - 1;
        end
        while aux_y >= 0 && y_filt(aux_y) - y_filt(aux_y - 1) >= 0
            aux_y = aux_y - 1;
        end
        b_on(i) = min(aux_x, aux_y);

        aux_x = b_off(i);
        aux_y = b_off(i);
        while aux_x <= L && x(aux_x) - x(aux_x + 1) >= 0
            aux_x = aux_x + 1;
        end
        while aux_y >= 0 && y(aux_y) - y(aux_y + 1) >= 0
            aux_y = aux_y - 1;
        end
        b_off(i) = max(aux_x, aux_y);
    end


    
    % Se o tempo entre 2 piscadas consecutivas for muito curto, fundo-as,
    % para não correr o risco de levar dados inúteis entre elas
    cnt = 1;
    while cnt > 0
    %     disp('while...')
        b = mov_props(b_on, b_off, data);
        itv_bl = length(b.interval);
        
        cnt = 0;
        for i = 1:itv_bl
            if b.interval(i) < par
                b_off(i)  = NaN;
                b_on(i+1) = NaN;
                cnt = cnt + 1;
            end
        end
        b_off(isnan(b_off)) = [];
        b_on(isnan(b_on))   = [];
    end

    b = mov_props(b_on, b_off, data);
end

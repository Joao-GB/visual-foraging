% As fixações são distinguidas das perseguições por uma medida de dispersão
function [fix, pur, rej_idx] = get_fixations_pursuits(data, nan_idx, disp_thr, m_f_len, v, v_thr)
% INPUTS:
    % data    : Posição do olho, em dva
    % fs      : frequência, em Hz
    % nan_idx : índices dos dados a serem desconsiderados para a busca de
    %            fixações/perseguições, mas que não foram removidos
    % disp_thr: limiar a partir do qual uma fixação passa a ser considerada
    %            pursuit, em dva
    % vel:    : vetor com as velocidades instantâneas, obtido via função
    %            get_saccades
    % vel_thr : limiar a partir do qual se considera que o vetor de velocidade
    %            descreve uma sacada, em dva/s
    %
    % OUTPUTS:
    % fixations : struct com info sobre fixações
    % pursuits  : struct com info sobre smooth pursuits
    L = length(data);
    aux = 1:L;
    aux(nan_idx) = NaN;
    s_idx = ~isnan(aux);
    % s de slow ou steady... Para confundir com sacada :P
    [s_on, s_off] = get_start_stop(s_idx);
    
    x = data(1,:);
    y = data(2,:);
    
    l = length(s_on);
    f_on = []; f_off = [];
    p_on = []; p_off = [];
    rej_idx = [];
    
    for i = 1:l
        len = s_off(i) - s_on(i);
%         disp(len)
        if len >= m_f_len || mean(v(s_on(i):s_off(i))) < v_thr
%             disp("Entrou!")
            x_range = max(x(s_on(i):s_off(i))) - min(x(s_on(i):s_off(i)));
            y_range = max(y(s_on(i):s_off(i))) - min(y(s_on(i):s_off(i)));
            dispersion = sqrt(x_range^2 + y_range^2);
            % Se a dispersão for baixa, temos uma fixação
            if dispersion <= disp_thr
                f_on  = [f_on  s_on(i)];
                f_off = [f_off s_off(i)];
            else
                p_on  = [p_on  s_on(i)];
                p_off = [p_off s_off(i)];
            end
        else
            rej_idx = [rej_idx s_on(i):s_off(i)];
        end
    end

    fix = mov_props(f_on, f_off, data);
    pur = mov_props(p_on, p_off, data);

end

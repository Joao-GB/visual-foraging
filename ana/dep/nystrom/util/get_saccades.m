% A primeira versão será baseada em Mitchell et al. (2014), com velocidade
% e aceleração calculadas ingenuamente
function [sac,vel,acc] = get_saccades(data, fs, nan_idx, vel_thr, acc_thr)
    % INPUTS:
    % data    : Posição do olho, em dva
    % fs      : frequência, em Hz
    % nan_idx : índices dos dados a serem desconsiderados para a busca de
    %            sacadas, mas que não foram removidos por sua possível
    %            contribuição espectral
    % vel_thr : limiar a partir do qual se considera que o vetor de velocidade
    %            descreve uma sacada, em dva/s
    % acc_thr : limiar a partir do qual se considera que o vetor de aceleração
    %            descreve uma sacada, em dva/s^2
    % method  : método de identificação de sacadas.
    %
    % OUTPUT:
    % saccades: struct como o blinks, cheio de informações relevantes sobre
    %           as sacadas

    lambda_on  = 5; % Onset da sacada é mais abrupto que offset
    lambda_off = 3;
    min_sac_dist = 50; % Se a diferença entre 2 sacadas for menor que 50 ms, funde-as
    min_sac_len  = 8;  % Se menor que uma microssacada, não passsa de ruído
    
    x = data(1,:);
    y = data(2,:);
    L = length(x);

    dt = 1/fs;      % O tempo entre 2 instantes é o inverso da frequência
    [vx, vy, ax, ay] = get_va1(x, y, fs, dt, nan_idx);

    % Fica tudo abaixo do limiar com essa coisa
%     [vx, vy, ax, ay] = get_va2(x, y, fs, nan_idx);

    vel = sqrt(vx.^2 + vy.^2); % Tamanho do vetor de velocidade
    acc = sqrt(ax.^2 + ay.^2); % Tamanho do vetor de aceleração

    
    % Os limiares são adaptativos, o máximo entre o dado pelo usuário e o
    % calculado com base no MAD dos dados
    acc_med = median(acc);
    acc_mad = median(abs(acc - acc_med));
    acc_thr_on = max(acc_thr, acc_med + lambda_on * acc_mad);
    acc_thr_off = max(acc_thr, acc_med + lambda_off * acc_mad);

    vel_med = median(vel);
    vel_mad = median(abs(vel - vel_med));
    vel_thr_on = max(vel_thr, vel_med + lambda_on * vel_mad);
    vel_thr_off = max(vel_thr, vel_med + lambda_off * vel_mad);


    % Como o limiar de offset é mais tolerante, primeiro identifico os
    % índices que estão acima dele
    v_aux = vel > vel_thr_off;
    a_aux = acc > acc_thr_off;
    aux = v_aux & a_aux;

    [sac_on, sac_off] = get_start_stop(aux);

    % Em seguida, adapto os onsets com base nos limiares de onset
    [sac_on, sac_off, not_sac] = change_onset(acc, vel, acc_thr_on, vel_thr_on, sac_on, sac_off, min_sac_dist, min_sac_len);

    sac = mov_props(sac_on, sac_off, data);
    % São índices a serem desconsiderados pois devem dizer respeito a
    % sacadas, mas não têm começo e/ou fim bem definidos devido ao ruído
    sac.idx_all = [sac.idx not_sac];
    sac.vel = vel;
    sac.vel_thr = vel_thr_off;

    % Desconsidero os intervalos que contêm NaN, pois pode ser que havia
    % alguma sacada nele
    l = size(sac.lims, 2);
    for i = 1:l-1
        if ~isempty(intersect(sac.lims(2, i):sac.lims(1, i+1), nan_idx))
            sac.interval(i) = NaN;
        end
    end
    sac.interval(isnan(sac.interval)) = [];

    sac.data = data;
    sac.data(:, setdiff(1:L, sac.idx)) = NaN;

end

% Posso achar a velocidade de dois jeitos
% 1) usando a versão ingênua de diff: get_saccade pode receber os dados
% com NaN que as bordas dos intervalos, quando subtraídas de NaN, viram NaN;
function [vx,vy,ax,ay] = get_va1(x, y, fs, dt, ni)
    fc = 50;
    order = 2;
    type = 'low';
    Wn = fc / (fs/2);
    [b, a] = butter(order,Wn,type);

    vx = diff(x) ./ dt;  % Velocidade horizontal
    vy = diff(y) ./ dt;  % Velocidade vertical

    % Filtro as velocidades antes de calcular a aceleração
    vx = filtfilt(b, a, vx);
    vy = filtfilt(b, a, vy);

    ax = diff(vx) ./ dt;
    ay = diff(vy) ./ dt;

    [vx, ax] = va_nans(vx, ax, ni, 0);
    [vy, ay] = va_nans(vy, ay, ni, 0);

end

% 2) Usando um filtro diferenciador: Get_saccade deve recebe data
    % inteiro e os instantes a serem desconsiderados em nan_idx. Após
    % aplicar o filtro (e obter o vetor de velocidades), transformo em NaN
    % os índices de nan_idx
    % Tenho que fazer uma função que aumenta os instantes de nan quando 
    % fazemos a diferença, para desconsiderar essas novas bordas no cáculo 
    % da velocidade e fazer a mesma coisa para obter acc_nan_idx:
    % z=zeros([1 length(x)]); z(blinks.idx) = NaN; z=diff(z); idx_v = find(isnan(z)); 
    % Verificação: length(idx_v)-length(blinks.idx)
function [vx,vy,ax,ay] = get_va2(x, y, fs, ni)
    
    % Define um filtro diferenciador como em Chen et al. (2021)
    fc = 54; N = 75; type = 'iv';
    b = diff_filter(fs, fc, N, type);

    vx = filter(b,1,x);
    vy = filter(b,1,y);
    ax = filter(b,1,vx);
    ay = filter(b,1,vy);
    
    [vx, ax] = va_nans(vx, ax, ni, 1);
    [vy, ay] = va_nans(vy, ay, ni, 1);
end


% Insere os NaNs no lugar certo para a velocidade e a aceleração
    % Preciso tornar NaN os índices das piscadas e seus vizinhos
    % antecessores (será que o filtro afeta a derivada dos posteriores?)
function [v1, a1] = va_nans(v, a, ni, method)
    if ~exist('method', 'var')
        method = 0;
    end
    v_ni = union(ni, ni-1);
    v_ni(v_ni <= 0) = [];   % Removo índices nulos (surgem se 1 \in ni)
    v(v_ni) = NaN;

    if method == 0
        % Calcular a diferença tira 1 elemento do vetor mas o alinhamento tem
        % poral importa. Solução: adiciono NaN no início do vetor, pois NaN < thr = 0,
        % logo vai realinhar com x e y sem afetar as comparações com limiar
        v1 = [NaN, v];
    
        a_ni = union(v_ni, v_ni-1);
        a_ni(a_ni <= 0) = [];   % Removo índices nulos (surgem se 1 \in ni)
        a(a_ni) = NaN;
        a1 = [NaN, NaN, a];
    else
        v1 = v;
        a1 = a;
    end
    
    % Para cada intervalo [start(i):stop(i)], encontro o primeiro instante
    % em que os limiares são superados e o redefino como novo onset
    % se o intervalo não superar esses limiares, é removido
end

function [start, stop, not_sac] = change_onset(a, v, at, vt, son, soff, min_d, min_l)
        not_sac = [];
        L = length(son);
        l = length(a);
        for i = 1:L
            % Intervalo que corresponde a uma sacada
            a_aux = a(son(i):soff(i));
            v_aux = v(son(i):soff(i));
            aux = a_aux > at & v_aux > vt;
            % Primeiro índice em que o limiar do onset é superado
            idx = find(aux, 1, 'first');
            % Se a sacada não superar esse limiar...
            if isempty(idx)
                son(i) = NaN;
                soff(i) = NaN;
            % ...ou se for antecedida ou seguida por NaN (vizinho de ruído),
            % a desconsidero
            elseif sum(isnan(a(son(i):max(son(i)-3, 1)))) > 0 || sum(isnan(a(soff(i):min(soff(i)+3, l)))) > 0
                not_sac = [not_sac son(i):soff(i)];
                son(i) = NaN;
                soff(i) = NaN;
            elseif soff(i) - son(i) < min_l
                son(i) = NaN;
                soff(i) = NaN;
            else
                son(i) = son(i) + idx;
            end
        end
        for i = 1:L-1
            if son(i+1) - soff(i) < min_d
                son(i+1) = son(i);
                son(i) = NaN;
                soff(i) = NaN;
            end
        end
        son(isnan(son)) = [];
        soff(isnan(soff)) = [];
        start = son;
        stop = soff;
    end
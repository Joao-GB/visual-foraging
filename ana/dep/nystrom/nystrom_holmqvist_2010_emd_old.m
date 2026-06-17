
% nystrom_get_psos não mais recalcula min_fix em função de fs...
function [eye_movs] = nystrom_holmqvist_2010_emd_old(data, fs, unit)
% No artigo, propõe primeiro filtrar, obtendo velocidade e aceleração, e
% depois mudar a unidade com um FATOR DE ESCALA... Completamente errado,
% mas vou fazer assim para obter resultados parecidos

% VERSÃO 2: primeiro mudo a escala do jeito certo, depois calculo a derivada
% com filtros SG


%% PARÂMETROS
% 0) Filtro para velocidade e aceleração
    a = 10;
    order = 2; 
    aux = ceil(2*a*fs/1000);
    len = aux + ~rem(aux,2); % Maior ou igual a 20 ms, mas deve ser ímpar

% 1) Mudança de unidade 
    pdr = 27.854;
    % dist = 45;

% 2) Piscadas
    min_itv = 0;    % Mínimo intervalo entre piscadas consecutivas, em ms
    pos_thr = 700;   % Atividade a partir da qual os dados saem da área de 
                     % interesse -- uso limiar conservador de 700 pixels
    vel_thr = 1000;  % Em dva/s
    acc_thr = 10000; % em dva/s^2

% 3) Sacadas
    min_sac = a;     % Duração mínima da sacada, em ms
    v_peak  = 120;   % Velocidade de pico inicial, deve ser entre 100 e 300 dva/s

% 4) Oscilações pós-sacádicas: sem parâmetros

% 5) Fixações
    min_fix = 40;    % Duração mínima da fixação, em ms


%% EXECUÇÃO
% 0) Filtro para velocidade e aceleração
    x = double(data(1,:));
    y = double(data(2,:));

   % A linha do meio de b é igual à 1a coluna de gSG (suavizar dados centrais)
    [~,gSG] = sgolay(order,len);
    gSG1 = gSG(:,2);
    gSG2 = gSG(:,3);

    vx = filtfilt(gSG1, 1, x);
    vy = filtfilt(gSG1, 1, y);

    ax = filtfilt(gSG2, 1, x);
    ay = filtfilt(gSG2, 1, y);

% 1) Mudança de unidade
    vel = fs*pdr*sqrt(vx.^2 + vy.^2); % Tamanho do vetor de velocidade
    acc = fs*pdr*sqrt(ax.^2 + ay.^2); % Tamanho do vetor de aceleração

% 2) Remoção de blinks
    % Uso critério de posição, aceleração e velocidade
    blinks = nystrom_get_blinks(data, vel, acc, min_itv, pos_thr, vel_thr, acc_thr);

    % Interpolo os dados de piscada com aproximação linear
    data(:, blinks.idx) = NaN;
    d = data;                                       % Rejeito as piscadas como NaNs...
    data(1,:) = fillmissing(data(1,:), 'linear');
    data(2,:) = fillmissing(data(2,:), 'linear');
    eye_movs.data = data;

% 3) Detecção de sacadas
    [saccades, thr_off, v_peak] = nystrom_get_saccades(data, vel, min_sac, v_peak, min_fix, fs);

% 4) Movimentos pós-sacádicos
    psos = nystrom_get_psos(data, vel, saccades, thr_off, v_peak, min_fix, fs);

% 5) Fixações
    nan_idx = [blinks.idx, saccades.idx, psos.idx];
    fixations = nystrom_get_fixations(data, nan_idx, min_fix, fs);



    % 6) Resultado final
    eye_movs.nan_data = d;
    eye_movs.vel = vel;
    eye_movs.acc = acc;
    eye_movs.blinks = blinks;
    eye_movs.fixations = fixations;
    eye_movs.saccades = saccades;
    eye_movs.psos = psos;

% 6) Fixações:
%     Exclusão se muito curtas, com o mesmo t_min
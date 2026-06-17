function [eye_movs] = eye_mov_detection(data, fs, unit)
% Exemplo:
% d541 = load("C:\Users\joaog\OneDrive\Área de Trabalho\USP\Mestrado\Projeto\data\free view\541\Anno_541_export.mat");
% data = d541.data.EYES;
% e = eye_mov_detection(data, 1000, 'pixel');

% Ordem dos processos:
% 0) Mudança de unidades
% 1) Filtração
% 2) Piscadas/ruído
% 3) Sacadas
% 4) Fixações/Smooth pursuit

%% PARÂMETROS a serem ajustados
% 0) Mudança de unidades
    pdr = 27.854;
    dist = 45;

% 1) Piscadas/ruído: deve ser feita antes da filtração, para não atrapalhar
    min_itv = 200;  % Mínimo intervalo entre piscadas (se grande, funde-as)
    bl_thr  = 30;   % Atividade a partir da qual os dados saem da área de 
                    % interesse (e.g., em pixels, a partir do momento que sair
                    % da tela)

% 2) Filtração
    fc = 40;        % Frequência de corte (50 Hz)
    window = 25;    % Janela do filtro de mediana
    order = 4;      % Para filtro Butterworth
    type = 'low';

% 3) Sacadas
    vel_thr  = 30;   % Em dva/s
    acc_thr  = 1000; % Em dva/s^2

% 4) Fixações/Smooth pursuit
    dispersion_thr = 1; % Em dva
    min_fix_len = 40;   % Comprimento mínimo para ser considerado algum 
                        %  movimento relevante


%% EXECUÇÃO

% 0) Mudança de unidades 
    if exist('unit', 'var') && strcmp(unit, "pixel")
        data = pixel_to_dva(data, 'dist', dist, 'pixel_dva_ratio', pdr);
    end


% 1) Piscadas/ruído
    % Identifico os momentos em que há piscadas e os salvo em blinks
    blinks = get_blinks(data, min_itv, bl_thr);

    % Interpolo os dados nos intervalos irrelevantes, apenas para não
    % afetar o desempenho dos filtros. Tentei pchip mas cria artefatos para
    % intervalos com bordas não-monotônicas
    data(:, blinks.idx) = NaN;
    d = data;                                       % Rejeito o as piscadas como NaNs...
    data(1,:) = fillmissing(data(1,:), 'linear');
    data(2,:) = fillmissing(data(2,:), 'linear');


% 2) Filtração
    data = eye_filter(data, fs, window, fc, order, type);
    d(~isnan(d)) = data(~isnan(d));


% 3) Sacadas
    saccades = get_saccades(data, fs, blinks.idx, vel_thr, acc_thr);


% 4) Fixações/Smooth pursuit
    nan_idx = [blinks.idx saccades.idx_all];
    min_fix_len = min_fix_len*1000/fs;
    [fixations, pursuits, rej_idx] = get_fixations_pursuits(data, nan_idx, dispersion_thr, min_fix_len, saccades.vel, vel_thr);
    
    d(:, rej_idx) = NaN;                            % ... E aqui rejeito o que não se enquadra em um dos tipos de movimento

    
% 5) Resultado final
    eye_movs.data = data;
    eye_movs.nan_data = d;
    eye_movs.blinks = blinks;
    eye_movs.fixations = fixations;
    eye_movs.pursuits = pursuits;
    eye_movs.saccades = saccades;
end
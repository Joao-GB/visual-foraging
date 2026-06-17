function [sac, ns_lims] = nystrom_get_saccades(data, fs, vel, nan_idx, options)
% NYSTROM_GET_SACCADES detecta os intervalos com sacadas em dados de
% movimento ocular conforme proposto em Nystrom & Holmqvist (2010)
% INPUTs:
% data    : Posição do olho, em dva
% vel     : velocidade dos olhos, em dva/s
% nan_idx : índices a serem ignorados para o cálculo de velocidade
% min_sac : duração mínima de uma sacada, em pts
% min_fix : duração mínima de uma fixação, em pts
% v_peak  : velocidade a ser ajustada iterativamente por nystrom_v_peak
% stop_thr: separação máxima entre limiares de pico de velocidade em
%           iterações consecutivas
% method  : 'mad' ou 'mean', para escolher como se calcula mu e sigma
%
% OUTPUTs:
% sac     : struct com informações sobre as sacadas

arguments
        data      {mustBeNumeric}
        fs        {mustBeNumeric}
        vel       {mustBeNumeric}
        nan_idx   {mustBeNumeric}
        options.minSac  {mustBeNumeric} = 10; 
        options.minFix  {mustBeNumeric} = 40;
        options.vPeak   {mustBeNumeric} = 100;
        options.method  string          = 'mean';
        options.stopThr {mustBeNumeric} = 1;
end
min_sac = options.minSac; min_fix = options.minFix; v_peak = options.vPeak; 
method  = options.method; stop_thr = options.stopThr;

% 1) Calcula o limiar de velocidade restrito aos índices de interesse
    vel_aux = vel; vel_aux(nan_idx) = NaN;
    v_peak = nystrom_v_peak2(vel_aux, v_peak, stop_thr, fs, min_sac, min_fix, method);
%     v_peak = nystrom_v_peak(vel_aux, v_peak, stop_thr, method);

% 2) Define preliminarmente os intervalos sacádicos
    sac = vel_aux >= v_peak;
    [s_on, s_off] = get_start_stop(sac);
    sac = mov_props(s_on, s_off);

% 3) Recalcula o início e o fim de cada sacada
    [sac, ns_lims] = nystrom_sac_on_n_off(data, fs, vel_aux, min_sac, min_fix, v_peak, sac, method);
end
function [eye_movs] = nystrom_holmqvist_2010_emd(eye_struct, options)
% NYSTROM_HOLMQVIST_2010_EMD implementa o algoritmo de detecção de 
% movimentos oculares proposto em Nystrom & Holmqvist (2010), que inclui 
% detecção de sacadas, fixações e oscilações pós-sacádicas.
%   NYSTROM_HOLMQVIST_2010_EMD(free_view_struct, "dataUnits", 'pixel',...
%                              "fieldData", 'EYES', "fieldFreq", 'Fs', "fieldPupil", 'PUPL',...
%                              "findEvents", {'blinks', 'saccades', 'fixations', 'glissades'},...
%                              "screenDist", 45, "screenRes", [1680 1050],
%                              "screenWidth", 47.376)
% OBS: assume que centro da tela corresponde ao (0,0)

arguments
    eye_struct struct
    options.dataUnits  string           = 'dva';
    options.fieldData  string           = 'EYES';
    options.fieldFreq  string           = 'Fs';
    options.fieldPupil string           = 'PUPL';
    options.findEvents cell             = {'blk', 'sac', 'fix', 'pso'};
    options.screenRes   {mustBeNumeric} = [1680 1050];
    options.screenWidth {mustBeNumeric} = 47.376;
    options.screenDist  {mustBeNumeric} = 45;
    options.pix2dvaDoneRight  {logical} = true;         % Opção para mudança de coordenadas do jeito correto
                                                        % ou do jeito proposto no artigo (fator de escala),
                                                        % que não é má aproximação para ângulos pequenos
    options.saveRawData {logical}       = false
    options.saveFiltData{logical}       = false
    options.saveMovProps{logical}       = false
end

    data    = eye_struct.(options.fieldData); pupil   = eye_struct.(options.fieldPupil);
    unit    = options.dataUnits;              fs      = eye_struct.(options.fieldFreq);
    ev_type = options.findEvents;             ev_type = change_ev_names(ev_type);
    all_ev_types = {'blinks', 'fixations', 'pursuits', 'saccades', 'psos'};

%% PARÂMETROS
% 0) Mudança de unidade
    width  = options.screenWidth(1);
    res    = options.screenRes(1);
    dist   = options.screenDist;

    exp_setup.frequency  = fs;
    exp_setup.resolution = options.screenRes;
    exp_setup.width      = options.screenWidth;
    exp_setup.distance   = dist;
    exp_setup.sessionID  = eye_struct.sessionID;


% 1) Piscadas/ruído: deve ser feita antes da filtração, para não atrapalhar
%    desempenho do filtro com os picos de ruído
    % (a) Filtro SG para vel e acc
    min_sac_ms   = 10;                         % Duração mínima de uma sacada, em ms

    l_SG_ms     = 2*min_sac_ms;                % Largura, em ms, da janela para filtro SG
    l_SG_pt     = ceil(l_SG_ms*fs/1000);       % Largura em pts 
    l_SG_pt     = l_SG_pt+ ~rem(l_SG_pt,2);
    o_SG        = 2;                           % Ordem do filtro SG

    % (b) Limiares que definem ruído
    vel_thr     = 1000;                        % Limiar de velocidade, em dva/s
    acc_thr     = 100000;                      % Limiar de aceleração, em dva/s^2

        % Os seguintes parâmetros não são dados no artigo
    pos_thr     = 30;                          % Limiar de posição, em dva
    pup_thr     = 0;                           % Limiar de tamanho de pupila
%     pos_thr     = pixel_to_dva(options.screenRes./2, "dist", dist, "res", res, "width", width);

    min_itv_ms  = 50;                           % Mínimo intervalo entre piscadas, em ms
    min_itv_pt  = ceil(min_itv_ms*fs/1000);    % em quantidade de pontos 


% 2) Sacadas e fixações
    min_sac_ms;                    %#ok<VUNUS> % Duração mínima de uma sacada, em ms
    min_sac_pt   = ceil(min_sac_ms*fs/1000);   % em quantidade de pontos

    min_fix_ms   = 40;                         % Duração mínima de uma fixação, em ms
    min_fix_pt   = ceil(min_fix_ms*fs/1000);

    v_peak       = 100;                        % Valor inicial para iterações feitas em nystrom_sac_on_n_off, em dva/s
    sep_peaks    = 1;                          % Valor para limitar o número de iterações, em dva/s
    method       = 'mean';                     % Método para calcular mu e sigma para limiares


%% EXECUÇÃO
%% 0) Mudança de unidades 
    if ~isa(data, 'double')
        data = double(data);
    end

    if strcmp(unit, "pixel")
        p2d_factor = pixel_to_dva(1, 'dist', dist, 'width', width, 'res', res);
        if options.pix2dvaDoneRight
            data = pixel_to_dva(data, 'dist', dist, 'width', width, 'res', res);
        else
            data = data.*p2d_factor;
        end
    end


%% 1) Piscada/ruído
    % (a) Obtém velocidade e aceleração com filtros SG
    [vel, acc] = nystrom_vel_n_acc(data, fs, o_SG, l_SG_pt);

    % (b) Identifica os instantes em que há piscadas e os salva em blinks
    blinks     = nystrom_get_blinks(data, pupil, vel, acc, "accThr", acc_thr, "velThr", vel_thr, "posThr", pos_thr, "pupThr", pup_thr, "itvThr", min_itv_pt);

    % (c) Usa instantes de piscadas como máscara
    data       = data_mask(data, blinks.idx, 'linear');

    % (d) Filtra os dados e atualiza vel e acc após linearização dos blinks
    %     em 1(c)
    [filt_data, vel, acc] = nystrom_sg_filt(data, fs, o_SG, l_SG_pt);


%% 2) Detecção de sacadas
    [saccades, non_sac] = nystrom_get_saccades(data, fs, vel, blinks.idx, "minSac", min_sac_pt, "minFix", min_fix_pt, "vPeak", v_peak, "method", method, "stopThr", sep_peaks);

    % Atualiza as piscadas com essas novas adições
    blk_plus_non_sac = merge_overlapping_events([non_sac blinks.lims]); blk_plus_non_sac = merge_adjacent_events(blk_plus_non_sac, max(1, min(min_fix_pt, min_sac_pt) - 1));
    [~, blk]  = subtract_overlapping_events(blk_plus_non_sac, blinks.lims);
    blinks = mov_props(blk);

%% 3) Movimentos pós-sacádicos
    [psos, saccades] = nystrom_get_psos(data, vel, saccades, min_fix_pt);

%% 4) Fixações
    blk_on  = blinks.lims(1,:); blk_on  = max(1, blk_on -1); 
    blk_off = blinks.lims(2,:); blk_off = min(length(data), blk_off+1);
    mod_blk_idx = mov_props([blk_on; blk_off]).idx;
    if any(strcmp(ev_type, all_ev_types{5}))
        nan_idx = sort([mod_blk_idx, saccades.idx, psos.idx]);
    else
        nan_idx = sort([mod_blk_idx, saccades.idx]);
    end

    [fix, unclassified] = nystrom_get_fixations(data, nan_idx, min_fix_pt);

%     unclassified = trim_overlapping_events(unclassified, [saccades.lims psos.lims]);
    blk_plus_unclass = reorder_events([blinks.lims, unclassified]); clear blinks;
    blk_plus_unclass = merge_adjacent_events(blk_plus_unclass, max(1, min(min_fix_pt, min_sac_pt) - 1));

    aux_fix = trim_overlapping_events(fix, blk_plus_unclass);
    fixations = mov_props(aux_fix);
    nFix = size(fixations.lims,2);
    pos = zeros(2,size(fixations.lims,2));
    for i=1:nFix
        pos(:,i) = mean(filt_data(:,fixations.lims(1,i):fixations.lims(2,i)), 2);
    end
    fixations.pos = pos;
    


%% 5) Resultado final
    blinks = mov_props(blk_plus_unclass);
    proper_lims = blinks.lims;                                             % Apesar de linealizar apenas os pontos de blinks.idx, para plotar
    proper_idx  = blinks.idx;                                              % preciso de sobreposição das bordas, evitando trajetórias isoladas.
                                                                           % Por isso esses fields a mais no struct blinks
    blk = proper_lims; 
    if ~isempty(blk)
        blk(1,:) = blk(1,:) - 1; if blk(1,1) <= 0; blk(1,1) = 1; end
        blk(2,:) = blk(2,:) + 1; if blk(2,end) > length(data); blk(2,end) = length(data); end
    end
    blinks = mov_props(blk); blinks.proper_lims = proper_lims; blinks.proper_idx = proper_idx;

    % Reaplica a máscara aos dados considerando a atualização das piscadas
    raw_data = data_mask(data, proper_idx); filt_data = data_mask(filt_data, proper_idx, 'linear'); 
%     nan_data = data_mask(filt_data, proper_idx);

    eye_movs.exp_setup     = exp_setup;
    eye_movs.emd_method    = 'Nystrom';
    if options.saveFiltData; eye_movs.data.filt_lin = filt_data; end
    if options.saveRawData;  eye_movs.data.raw_nan  = raw_data;  end 
    eye_movs.data.unit     = 'dva';
    eye_movs.data.pos_lims = pos_thr;
    eye_movs.time_unit     = 'ms';
    eye_movs.vel           = vel;
    eye_movs.acc           = acc;

    if ~options.saveMovProps
        fields_to_remove = {'duration', 'idx', 'interval', 'proper_idx'};
        blinks = rmfield(blinks, fields_to_remove); fields_to_remove(end) = [];
        saccades = rmfield(saccades, fields_to_remove);
        fixations = rmfield(fixations, fields_to_remove);
        psos = rmfield(psos, fields_to_remove);
    end

    if any(strcmp(ev_type, all_ev_types{1})); eye_movs.blinks = blinks;       end
    if any(strcmp(ev_type, all_ev_types{2})); eye_movs.fixations = fixations; end
    if any(strcmp(ev_type, all_ev_types{4})); eye_movs.saccades = saccades;   end
    if any(strcmp(ev_type, all_ev_types{5})); eye_movs.psos = psos;           end

    aux = []; idx_aux = [];
    em_ev = fieldnames(eye_movs);
    for i = 1: length(all_ev_types)
        if any(strcmp(em_ev, all_ev_types{i}))
            aux = [aux eye_movs.(all_ev_types{i}).lims]; %#ok<*AGROW> 
            idx_aux = [idx_aux ones(1, size(eye_movs.(all_ev_types{i}).lims,2))*(i-1)];
        end
    end
    eye_movs.all_lims = [aux; idx_aux];
    eye_movs.all_lims = reorder_events(eye_movs.all_lims);
end
function eye_movs = microsacc_emd(eye_struct, options)
% OBS: assume que centro da tela corresponde ao (0,0)
arguments
    eye_struct struct
    options.dataUnits  string           = 'pixel'
    options.fieldData  string           = 'EYES'
    options.fieldFreq  string           = 'Fs'
    options.fieldPupil string           = 'PUPL'
    options.findEvents cell             = {'blk', 'sac', 'fix'}
    options.screenRes   {mustBeNumeric} = [1680 1050];
    options.screenWidth {mustBeNumeric} = 47.376;
    options.screenDist  {mustBeNumeric} = 45;
    options.saveRawData {logical}       = false
    options.saveFiltData{logical}       = false
    options.saveMovProps{logical}       = false
end

data    = eye_struct.(options.fieldData);
fs      = eye_struct.(options.fieldFreq);
pupil   = eye_struct.(options.fieldPupil);
unit    = options.dataUnits;
ev_type = options.findEvents;
ev_type = change_ev_names(ev_type);
all_ev_types = {'blinks', 'fixations', 'pursuits', 'saccades', 'psos'};

%% PARÂMETROS
% 0) Mudança de unidades
    width  = options.screenWidth(1);
    res    = options.screenRes(1);
    dist   = options.screenDist;

    exp_setup.frequency  = fs;
    exp_setup.resolution = options.screenRes;
    exp_setup.distance   = dist;


% 1) Piscadas/ruído: deve ser feita antes da filtração, para não atrapalhar
%    desempenho do filtro com os picos de ruído
    % (a) Filtro SG para vel e acc
    l_SG_ms     = 20;                          % Largura, em ms, da janela para filtro SG
    l_SG_pt     = ceil(l_SG_ms*fs/1000);       % Largura em pts 
    l_SG_pt     = l_SG_pt+ ~rem(l_SG_pt,2);
    o_SG        = 2;                           % Ordem do filtro SG

    % (b) Limiares que definem ruído
    vel_thr     = 20000;                       % Limiar de velocidade, em dva/s
    acc_thr     = 200000;                      % Limiar de aceleração, em dva/s^2
    pos_thr      = 30;                         % Limiar de posição, em dva
    pup_thr     = 0;                           % Limiar de tamanho de pupila

    min_itv_ms  = 50;                          % Mínimo intervalo entre piscadas, em ms
    min_itv_pt  = ceil(min_itv_ms*fs/1000);    % em quantidade de pontos


% 2) Sacadas e fixações
    % Parâmetros do próprio algoritmo microsacc
    MINDUR = 5;        % Minimum duration (number of samples)
    VTHRES = 4;        % Velocity threshold
    min_sac_ms   = 10;                         % Duração mínima de uma sacada, em ms
    min_sac_pt   = ceil(min_sac_ms*fs/1000);   % em quantidade de pontos 
%     min_sac_ms_2 = 5;                          % Outra duração mínima de sacada, em ms (cf. mod_CluterFix)

    min_fix_ms   = 40;                         % Duração mínima de uma fixação, em ms
    min_fix_pt   = ceil(min_fix_ms*fs/1000);


%% EXECUÇÃO
if ~isa(data, 'double')
    data = double(data);
end

%% 0) Mudança de unidades 
    if strcmp(unit, "pixel")
        data = pixel_to_dva(data, 'dist', dist, 'width', width, 'res', res);
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
    [filt_data, vel, acc] = nystrom_sg_filt(data, fs, o_SG, l_SG_pt, 1);
    vx = vel.vx; vy = vel.vy; vel = vel.vel; acc = acc.acc;


%% 2) Sacadas e Fixações
% Caso eu quisesse trabalhar com a velocidade do próprio algoritmo
    aux_vel = vecvel(data',fs,2);
    [sac, ~] = microsacc(data',aux_vel, VTHRES, MINDUR);
%     [sac, ~] = microsacc(data',[vx' vy'], VTHRES, MINDUR);

    sac_lims = sac(:, 1:2)';
    sac_lims = merge_adjacent_events(sac_lims, ceil(min_fix_pt/2));
    [sac, non_sac] = subtract_overlapping_events(sac_lims, blinks.lims);
    sac = mov_props(sac);

    % (c) Atualiza as piscadas com essas novas adições
    blk_plus_non_sac = merge_overlapping_events([non_sac blinks.lims]); blk_plus_non_sac = merge_adjacent_events(blk_plus_non_sac, max(1, min(min_fix_pt, min_sac_pt) - 1));
    [~, blk]  = subtract_overlapping_events(blk_plus_non_sac, blinks.lims);
    blinks = mov_props(blk);

    if ~isempty(blinks.lims)
        blk_on  = blinks.lims(1,:); blk_on  = max(1, blk_on -1); 
        blk_off = blinks.lims(2,:); blk_off = min(length(data), blk_off+1);
        mod_blk_idx = mov_props([blk_on; blk_off]).idx;
    else
        mod_blk_idx = [];
    end

    nan_idx = sort([mod_blk_idx, sac.idx]);

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
    saccades = sac; 
    % Recalcula a velocidade para que seja dada pelos métodos dele
    vel = sqrt(aux_vel(:,1).^2+aux_vel(:,2).^2); vel = vel';
    [v_peaks, amp] = get_v_peak_n_amp(data, vel, saccades.lims, 'max'); 
    saccades.v_peaks = v_peaks.v_peak; saccades.v_peak_idx = v_peaks.idx;  saccades.amplitude = amp;

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
    data = data_mask(data, proper_idx, 'linear'); raw_data = data_mask(data, proper_idx); 
    filt_data = data_mask(filt_data, proper_idx, 'linear'); 
%     nan_data = data_mask(filt_data, proper_idx);

    eye_movs.exp_setup     = exp_setup;
    eye_movs.emd_method    = 'Microsacc';
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
    end
        
    if any(strcmp(ev_type, all_ev_types{1})); eye_movs.blinks = blinks;       end
    if any(strcmp(ev_type, all_ev_types{2})); eye_movs.fixations = fixations; end
    if any(strcmp(ev_type, all_ev_types{3})); eye_movs.pursuits = pursuits;   end
    if any(strcmp(ev_type, all_ev_types{4})); eye_movs.saccades = saccades;   end
    if any(strcmp(ev_type, all_ev_types{5})); eye_movs.psos = pso;            end

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
    % Para teste, sum(eye_movs.all_lims(2,:)-eye_movs.all_lims(1,:)) deve
    % dar L - 1, em que L = size(data, 2) (soma telescópica)
end
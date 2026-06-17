function [pso, sac] = nystrom_get_psos(data, vel, sac, min_fix)
    sac_lims = sac.lims; sac_amp  = sac.amplitude; thr_off  = sac.thr_off; 
    sac_vp_thr = sac.vp_thr; sac_vpeaks = sac.v_peaks; sac_vp_idx = sac.v_peak_idx;
    sac.thr_off = thr_off;
    L = length(data); l = length(sac_lims(2,:));

    sac_peaks_idx = zeros(1, l);
    for i =1:l
        [~,aux] = findpeaks(vel(sac_lims(1,i):sac_lims(2,i)), "NPeaks",1);
        sac_peaks_idx(i) = sac_lims(1,i) - 1 + aux;
    end

% 1) Após cada sacada, procura numa janela de tamanho min_fix atividade
% oscilatória
    p_lims = nan(2, l); %     type = zeros(1, l);
    i = 1;
    while i <= l
%         if i >= 25
%             pause(.5);
%         end
        step = 0;
        post_s_itv = sac_lims(2,i):min(sac_lims(2,i) + min_fix - 1,L);

% 2) Verifica se é pso de baixa velocidade
    % (a) Obtém os os índices em que a atividade pós-sacádica é rápida
        p_bol = vel(post_s_itv) >= thr_off(i);
        [~,p_low_off] = get_start_stop(p_bol);
    % (b) Os índices de alta velocidade não podem estar concentrados perto
    %     do fim de pos_s_itv
        p_low_off(p_low_off == length(post_s_itv)) = [];

    % (c) Se houver algum intervalo rápido em post_s_itv, define o pso
        if ~isempty(p_low_off)
            p_lims(1,i) = sac_lims(2,i);
            p_lims(2,i) = p_lims(1,i) + p_low_off(end);
        end

% 3) Verifica se é pso de alta velocidade (era sacada)
        if i+1 <= l
    % (a) Verifica se alguma sacada seguinte à atual está próxima (i.e.,
    %     com o pico dentro do intervalo pos_s_itv
            p_high_off = find(sac_peaks_idx((i+1):end) <= post_s_itv(end));

            if any(p_high_off)
    % (b) Se tal sacada seguinte tiver amplitude 2x menor que a da sacada 
    % atual, é pso de alta velocidade
                idx = p_high_off(end);
                if sac_amp(i) >= 2*max(sac_amp((i+1):(i+idx))) && ((~isnan(p_lims(2,i)) && sac_lims(2,i+idx)>p_lims(2,i)) || isnan(p_lims(2,i)))
                    fprintf('Exclusão: intervalo de PSO')
                    for k = 1:idx
                        fprintf(', %d a %d', sac_lims(1,i+k), sac_lims(2,i+k));
                    end
                        fprintf('.\n')
                    p_lims(1,i) = sac_lims(2,i);
                    p_lims(2,i) = sac_lims(2,i+idx);
                    sac_lims(1,i+(1:idx)) = NaN; sac_lims(2,i+(1:idx)) = NaN;
                    sac_amp(i+(1:idx)) = NaN; thr_off(i+(1:idx)) = NaN; sac_vpeaks(i+(1:idx)) = NaN; sac_vp_idx(i+(1:idx)) = NaN;
                    step = idx;
                else
                end
            end
        end
    
    % (c) Determina o fim do pso, caso tenha sido identificado algum
        if ~isnan(p_lims(2,i))
            while p_lims(2,i) < length(vel) && (vel(p_lims(2,i)) > thr_off(i) || diff(vel(p_lims(2,i)+(0:1))) < 0)
                    p_lims(2,i) = p_lims(2,i)+1;
            end
        end
        i = i + step + 1;
    end
% 4) Define os structs de psos
    p_on = p_lims(1,:); p_off = p_lims(2,:);
    p_on(isnan(p_on))     = [];
    p_off(isnan(p_off))   = []; %     type(isnan(type))   = [];

    pso = mov_props([p_on; p_off]); %     pso.type = type;

% 5) Remove as sacadas que viraram psos
    ns_lims_on = sac.lims(1,isnan(sac_lims(1,:))); ns_lims_off = sac.lims(2,isnan(sac_lims(2,:)));
%     ns_lims = [ns_lims_on; ns_lims_off];
    sac = mov_props(sac_lims);
    sac_amp(isnan(sac_amp)) = []; thr_off(isnan(thr_off)) = []; sac_vpeaks(isnan(sac_vpeaks)) = []; sac_vp_idx(isnan(sac_vp_idx)) = [];
    
    sac.amplitude = sac_amp; sac.thr_off = thr_off; 
    sac.vp_thr = sac_vp_thr; sac.v_peaks = sac_vpeaks; sac.v_peak_idx = sac_vp_idx;
end
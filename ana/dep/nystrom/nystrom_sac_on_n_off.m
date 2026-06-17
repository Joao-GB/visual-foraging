function [sac, non_sac_lims] = nystrom_sac_on_n_off(data, fs, vel_nan, min_sac, min_fix, vp_thr, sac, method, sep_thr)
    % Parâmetros de peso para limiar de offset
    og_sac = sac;
    a = .7; b = 1 - a;
    half_sac = ceil(min_sac/2);
    
    s_on = sac.lims(1,:); s_off = sac.lims(2,:);
    L = length(vel_nan);


% 1) Se sep_thr é dado como entrada, vp_thr não foi calculado ainda.
%        (se estiver usando via nystrom_get_saccades, basta não passar
%         sep_thr que esse reajuste não será feito)
    if exist('sep_thr', 'var')
        vp_thr = nystrom_v_peak2(vel_nan, vp_thr, sep_thr, fs, min_sac, min_fix, method);
    %     vp_thr  = nystrom_v_peak(vel, vp_thr, sep_thr);
    end


% 2) Obtém o pico de velocidade de cada potencial sacada
    v_peak     = nan(1, length(s_on));
    v_peak_idx = nan(1, length(s_on));
    for i=1:length(s_on)

        try
            % (a) Procura por picos no intervalo
            [s_vp, si_vp] = findpeaks(vel_nan(s_on(i):s_off(i)), SortStr='descend', NPeaks=1);
            v_peak(i)     = s_vp;
            v_peak_idx(i) = si_vp + s_on(i) - 1;
            [~, si_vM] = max(vel_nan(s_on(i):s_off(i)));
            
            % (b) Se o pico não coincidir com o máximo, não é sacada pura
            if si_vp ~= si_vM
                disp(['Exclusão: Pico diferente de máximo, ' num2str(s_on(i)) ' a ' num2str(s_off(i))])
                v_peak(i) = NaN; v_peak_idx(i) = NaN; 
                s_on(i) = NaN; s_off(i) = NaN;
            end
        catch ME
            % (c) Se intervalo muito curto ou sem picos, não deve ser sacada
            if ismember(ME.identifier, {'signal:findpeaks:emptyDataSet', 'MATLAB:matrix:singleSubscriptNumelMismatch'})
                disp(['Exclusão: intervalo sem pico, ' num2str(s_on(i)) ' a ' num2str(s_off(i))])
                s_on(i) = NaN; s_off(i) = NaN;
            else; rethrow(ME);
            end
        end
    
    
% 3) Se sep_thr é dado como entrada, encurta as sacadas antes de procurar
%    seu início e fim
        if exist('sep_thr', 'var') && ~isnan(s_on(i)) && s_off(i) - s_on(i) > min_sac
            s_off(i) = s_on(i) + (si_vp - 1) + 1;
            s_on(i)  = s_on(i) + (si_vp - 1) - 1;
        end
    end


% 4) Obtém o limiar (global) de início de sacada
    [mu, sigma] = nystrom_mu_sigma_method(vel_nan(vel_nan < vp_thr), method);
    thr_on = mu + 3*sigma;


% 5) Obtém o início, o fim e a amplitude de cada sacada
    amp = nan(1, length(s_on)); thr_off = zeros(1, length(s_on));

    s_idx = find(~isnan(s_on));
    for i = s_idx

    % (a) Recua o onset até que as duas condições de interesse sejam
    %     satisfeitas. No entanto, para se chegar num NaN
        while s_on(i) - 1 >= 1 && ~isnan(vel_nan(s_on(i) - 1)) && ~(vel_nan(s_on(i)) < thr_on && (vel_nan(s_on(i) - 1) - vel_nan(s_on(i))) >= 0)
            s_on(i) = s_on(i) - 1;
        end
    
    % (b) Encontra os min_fix pontos não nulos antes do início da sacada, a
    %     saber, que não são sacadas nem NaNs
        sac_bol = zeros(1, L); sac_bol(sac.idx) = 1;
        pre_s_itv = find(~isnan(vel_nan(1:s_on(i)-1)) & ~sac_bol(1:s_on(i)-1), min_fix,'last');

        
    % (c) Obtém o limiar (local) de fim de sacada
        [mu_t, sigma_t] = nystrom_mu_sigma_method(vel_nan(pre_s_itv), method);
        thr_t = mu_t + 3*sigma_t;

        % Se o limiar local for muito alto, não é usado em thr_off
        if ~isnan(thr_t) && thr_t < vp_thr;  thr_off(i) = a*thr_on + b*thr_t;
        else;                               thr_off(i) = thr_on;
        end

    % (d) Avança o offset até que as duas condições de interesse sejam
    %     satisfeitas. No entanto, para se chegar num NaN
        while s_off(i) + 1 <= L && ~isnan(vel_nan(s_off(i) + 1)) && ~(vel_nan(s_off(i)) <= thr_off(i) && (vel_nan(s_off(i)) - vel_nan(s_off(i) + 1)) <= 0)
            s_off(i) = s_off(i) + 1;
        end

    % (e) Exclui a potencial sacada se:
        %   i. se não houver período pré-sacádico longo o suficiente
        %      (i.e., apenas se a sacada estiver muito perto do início)
        %  ii. se o intervalo anterior à candidata sacada for de 
        %      alta velocidade;
        % iii. se a sacada for muito curta.
        if length(pre_s_itv) <=1 || mu_t > vp_thr ||  s_off(i) - s_on(i) < min_sac
            disp(['Exclusão: intervalo curto ou com baixa velocidade, ' num2str(s_on(i)) ' a ' num2str(s_off(i))])
            s_on(i)     = NaN; s_off(i)    = NaN;
            thr_off(i)  = NaN; v_peak(i)   = NaN;  v_peak_idx(i) = NaN; 
        end

     % (f) Se duas sacadas se sobrepuserem, uma delas é excluída
        if i > 1 && s_on(i) < s_off(i-1)
            s_on(i)    = min(s_on(i), s_on(i-1));       s_off(i)   = max(s_off(i), s_off(i-1));
            thr_off(i) = max([thr_off(i), thr_off(i-1)]); [v_peak(i), idx]  = max([v_peak(i-1), v_peak(i)]); 
            v_peak_idx(i) = v_peak_idx(i + idx - 2);
            s_on(i-1)     = NaN; s_off(i-1)    = NaN;
            thr_off(i-1)  = NaN; v_peak(i-1)   = NaN; v_peak_idx(i-1) = NaN; 
            amp(i-1)      = NaN;
        end
        amp(i) = get_distance(data, s_on(i), s_off(i));
    % (g) Obtém o struct da sacada com os novos limites calculados
    %     (para calcular sac_bol, que usa sac.idx)
        sac    = mov_props(s_on, s_off);
    end

    non_sac_lims = [og_sac.lims(1, isnan(s_on)); og_sac.lims(2, isnan(s_off))];
    nan_idx = find(isnan(s_on));
    s_on(nan_idx)    = [];
    s_off(nan_idx)   = [];
    thr_off(nan_idx) = [];
    v_peak(nan_idx)  = []; v_peak_idx(nan_idx)  = []; 
    amp(nan_idx)     = [];

    sac = mov_props([s_on; s_off]);
    sac.vp_thr  = vp_thr; sac.v_peaks = v_peak; sac.v_peak_idx = v_peak_idx; sac.amplitude = amp;
    sac.thr_off = thr_off;
    
    function d = get_distance(dat, i1, i2)
        if isnan(i1)
            d = NaN;
        else
            d = sqrt( (dat(1,i1)-dat(1,i2))^2 + (dat(2,i1)-dat(2,i2))^2 );
        end
    end
end
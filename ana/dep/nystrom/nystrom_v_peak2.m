function v_peak = nystrom_v_peak2(vel, v_0, stop, fs, min_sac, min_fix, method)
% Calcula o limiar de velocidade a ser utilizado para detecção de sacadas,
% modificado de Nystrom & Holmqvist (2010), a partir da função doOptimize,
% cf. https://github.com/dcnieho/NystromHolmqvist2010/blob/master/function_library/estimateThresholds.m
% vel   : vetor de velocidades
% v_0   : velocidade-chute inicial a ser ajustada iterativamente. Se -1,
%         usa um valor incial também data-driven
% stop  : condição de parada da iteração, em dva/s
% method: 'median' ou 'mean'

shorten_fix   = ceil(min_sac/3000 * fs);
if v_0 == -1
    v_aux = rmoutliers(vel, 'median'); v_0 = median(v_aux) + 3*std(v_aux);
end

count = 0;
time_stop = 5000;

v_peak_prev = Inf;
v_peak = v_0;
while abs(v_peak - v_peak_prev) >= stop && count < time_stop
    idx_v_below = vel <= v_peak;
    [vb_on, vb_off] = get_start_stop(idx_v_below);

    idx_vb_long = vb_off-vb_on >= min_fix;
    vb_on    = vb_on(idx_vb_long);
    vb_off   = vb_off(idx_vb_long);
    
    % shrink them as done in Nyström's version, to make sure we don't
    % catch the data that is still during the saccade
    vb_on    = vb_on  + shorten_fix;
    vb_off   = vb_off - shorten_fix;

    aux = mov_props(vb_on, vb_off);

    
    v_peak_prev = v_peak;
    [mu, sigma] = nystrom_mu_sigma_method(vel(aux.idx), method);
    v_peak = mu + 6*sigma;

    count = count + 1;
end
end
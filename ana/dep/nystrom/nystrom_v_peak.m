function v_peak = nystrom_v_peak(vel, v_0, stop, method)
% Calcula o limiar de velocidade a ser utilizado para detecção de sacadas,
% conforme Nystrom & Holmqvist (2010).
% vel   : vetor de velocidades
% v_0   : velocidade-chute inicial a ser ajustada iterativamente. Se -1,
%         usa um valor incial também data-driven
% stop  : condição de parada da iteração, em dva/s
% method: 'median' ou 'mean'
if nargin < 4; method = 'median'; end

if v_0 == -1
    v_aux = rmoutliers(vel, 'median'); v_0 = mean(v_aux) + 3*std(v_aux);
end


count = 0;
time_stop = 5000;

v_peak_prev = Inf;
v_peak = v_0;
while abs(v_peak - v_peak_prev) >= stop && count < time_stop
    v_below = vel(vel <= v_peak);
    [mu, sigma] = nystrom_mu_sigma_method(v_below, method);
    v_peak_prev = v_peak;
    v_peak = mu + 6*sigma;

    count = count + 1;
%     disp(v_peak)
end
% disp(count)
end
function trials_heteroscedasticity(data, ev, name, rad_or_ang)
% data      : dados da posição do olho
% ev        : struct com informações sobre o evento (obtido de mix_emd, por
%             exemplo
% name      : nome do evento analisado, para aparecer no título do gráfico
% rad_or_ang: nome do modo a ser utilizado
    i_on  = ev.lims(1,:);
    i_off = ev.lims(2,:);

% Dependências: circ_stats
    x = data(1,:);
    y = data(2,:);

    r = sqrt(x.^2+y.^2);
    a = atan(y./x);

    l = length(i_on);
    r_means = zeros([1, l]); r_stds  = zeros([1, l]);

    a_means = zeros([1, l]); a_stds  = zeros([1, l]);

    % Para cada segmento de dado (e.g. para cada fixação), obtém a distância
    % média (e sua variância) e o ângulo
    for i = 1:l
        r_seg      = r(i_on(i):i_off(i));
        r_means(i) = mean(r_seg);
        r_stds(i)  = std(r_seg);

        a_seg      = a(i_on(i):i_off(i));
        a_means(i) = rad2deg(circ_mean(a_seg'));
        a_stds(i)  = rad2deg(circ_std(a_seg'));
    end
    [r_means, r_stds] = is_valid(r_means, r_stds);
    [a_means, a_stds] = is_valid(a_means, a_stds);

%     rad = {'rad', 'radius', 'radii', 'r'};
%     ang = {'ang', 'angle', 'angles', 'a'};

    figure; 
    subplot(2,1,1)
    % Plota apenas a variância em função da distância como um scatter plot
    scatter(r_means, r_stds)

    % Tenta fazer regressão linear da variância das sacadas em função da
    % distância média 
    P = lin_reg_fit(r_means, r_stds);
    yfit = polyval(P,r_means);
    hold on; plot(r_means,yfit,'r-');
        
    xlabel("Mean radial eye position (dva)"); ylabel("Standard deviation");
    title(['Radial eye variability during ' name])

    eqn = string(" Linear: y = " + P(1)) + "x + " + string(P(2));
    text(min(r_means),max(r_stds),eqn,"HorizontalAlignment","left","VerticalAlignment","top")
        

    subplot(2,1,2);
    % Plota apenas a variância em função da da posição angular como um scatter plot  
    % Removo outliers do y (a_stds)
    out = rem_out(a_stds);
    a_means = a_means(~out);
    a_stds = a_stds(~out);
    
    polarscatter(a_means, a_stds);
    title(['Angular eye variability during ' name])

end

function p_final = lin_reg_fit(x, y)

% Fit/ajuste inicial
p_init = polyfit(x, y, 1);
y_pred = polyval(p_init, x);
residuals = y - y_pred;

out = rem_out(residuals);

% Ajuste final
x1 = x(~out);
y1 = y(~out);
p_final = polyfit(x1, y1, 1);
end

function out = rem_out(y)
    lambda = 3;
    center = median(y); % Robust center estimate
    abs_dev = abs(y- center);
    mad = median(abs_dev); % Median Absolute Deviation
        
    % Define threshold (3 MADs is common for outliers)
    threshold = lambda * mad;
    out = abs_dev > threshold;
end

function [x1, y1] = is_valid(x, y)
    valid = ~(isnan(x) | isnan(y) | isinf(x) | isinf(y));
    x1 = x(valid);
    y1= y(valid);
end

% OUTRAS TENTATIVAS...
    % Tentativa 1: regressão linear com polyfit
    % P = polyfit(r_means,r_stds,1);
    % Tentativa 2: regressão linear sem outliers
%         P = robustfit(r_means,r_stds);
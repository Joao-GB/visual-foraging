function [f_data, vel, acc] = nystrom_sg_filt(data, fs, order, len, verbose)
% NYSTROM_SG_FILT obtém os vetores de velocidade, aceleração e dados 
% filtrados para dados de posição ocular utilizando filtro de Savitzky-Golay 
% conforme Nystrom & Holmqvist (2010)
%   NYSTROM_SG_FILT(eye_data, 1000, 2, 21)
% 
% INPUTS:
%   data   : vetor 2 x L com dados de posição do olho, em dva%   event_type  : string, e.g. 'saccade'
%   fs     : frequência de amostragem
%   order  : ordem do filtro SG
%   len    : comprimento do filtro SG
%   verbose: valor para mudar o tipo de output da função (cf. abaixo)
% 
% OUTPUT:
%   f_data : dados filtrados
%   vel    : array com velocidade (se verbose não for dado) ou struct com
%            vel, vx e vy (se verbose for dado como input)
%   acc    : array com aceleração (se verbose não for dado) ou struct com
%            acc, ax e ay (se verbose for dado como input)

    x = data(1,:);
    y = data(2,:);

    pad_len = ceil(len/2);
    x_pad = apply_padding(x, pad_len, 'add', 'symmetric');
    y_pad = apply_padding(y, pad_len, 'add', 'symmetric');

    % Filtros de velocidade e aceleração
    % A linha do meio de b é igual à 1a coluna de gSG (suavizar dados centrais)
    [~,gSG] = sgolay(order,len);

    % Preciso inverter o array de coeficientes da transfer function devido
    % à maneira como a convolução é implementada. Se usasse o filtro para
    % obter o produto interno entre gSG1 e o vetor de tempo lagged, não
    % precisaria do flipud.Para 'len' ímpar (recomendado), flipud equivale
    % a inverter o sinal, já que gSG1 é antissimétrico
    dt = 1/fs;
    gSG0 = flipud(gSG(:,1));     % coeficientes para suavizar os dados
    gSG1 = flipud(gSG(:,2)) / dt;     % coeficientes de 1a derivada
    gSG2 = flipud(gSG(:,3)) / dt^2;   % coeficientes de 2a derivada
    
    x      = conv(x_pad, gSG0, 'same'); x = apply_padding(x, pad_len, 'remove');
    y      = conv(y_pad, gSG0, 'same'); y = apply_padding(y, pad_len, 'remove');
    f_data = [x; y];

    vx = conv(x_pad, gSG1, 'same');
    vy = conv(y_pad, gSG1, 'same');

    ax = conv(x_pad, gSG2, 'same');
    ay = conv(y_pad, gSG2, 'same');

    vel = sqrt(vx.^2 + vy.^2); % Velocidade absoluta
    vel = apply_padding(vel, pad_len, 'remove');
    acc = sqrt(ax.^2 + ay.^2); % Aceleração absoluta
    acc = apply_padding(acc, pad_len, 'remove');

    if nargin >=5
       vx = apply_padding(vx, pad_len, 'remove');
       vy = apply_padding(vy, pad_len, 'remove');
       ax = apply_padding(ax, pad_len, 'remove');
       ay = apply_padding(ay, pad_len, 'remove');
       v.vel = vel;
       v.vx  = vx;
       v.vy  = vy;
       a.acc = acc;
       a.ax  = ax;
       a.ay  = ay;
       vel = v;
       acc = a;
    end
end
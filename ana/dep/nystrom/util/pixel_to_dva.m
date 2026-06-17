function d = pixel_to_dva(data, options)
% Converte data em pixels para d em dva (degrees of visual angle), usando 2
% possíveis entradas distintas
% 1) Se for dado pixel_dva_ratio, i.e., a quantidade de pixels que há em
%    exatamente 1 dva;
% 2) Se forem passadas a largura da tela (width) e sua resolução (res)
% INPUT
% data:            dados de eye-tracking em pixels;
% dist:            distância entre o monitor e o sujeito, em cm. Padrão: 45;
% pixel_dva_ratio: quantos pixels são compreendidos por 1° (lembre: a relação
%                   não é linear, então para 2° não será o dobro);
% width:           largura da tela, em cm;
% res:             resolução da tela, em pixels (e.g. [1920, 1080] ou
%                   [1680, 1050], mas só importa a largura).

% OUTPUT:
% d:               dados em dva.

% EXEMPLO: d2 = pixel_to_dva(data, 'dist', 45, 'width', 47.376, 'res', 1680);
%          d3 = pixel_to_dva(data, 'dist', 45, 'pixel_dva_ratio', 27.854);
    
%     options =  dva_n_pixel_set_options(varargin);
    arguments
        data;
        options.dist (1,1) double = 45
        options.pixel_dva_ratio (1,1) double = NaN
        options.width (1,1) double = NaN
        options.res (1,1) double = NaN
    end
    
    if ~isnan(options.pixel_dva_ratio)
        pdr = options.pixel_dva_ratio;
    
        % Comprimento, em cm, da região da tela coberta por 1 dva
        % Tangente de 1 dva é cateto oposto (a ser descoberto) sobre
        % adjacente (distância)
        l = options.dist*tan(pi/180);
        cpr = l/pdr; % cpr = cm_pixel_ratio (i.e., largura de 1 pixel em cm)
    
        d = pixel2rad2dva(data, cpr, options.dist);
    
    elseif ~isnan(options.width) && ~isnan(options.res(1))
    
        cpr = options.width/options.res(1); % cpr = cm_pixel_ratio 
        d = pixel2rad2dva(data, cpr, options.dist);
    
    else
        disp("Passe algum argumento adicional para transformar seus dados.")
        d = data;
    end

end

% Converte pix em radiano e depois em dva
function deg = pixel2rad2dva(pix, c, dist)
    rad = atan(pix.*(c/dist)); % Radiano é arctan de cateto oposto/adjacente;
                               % Dados em pixel são multiplicados pelo
                               % tamanho de cada pixel, obtendo-se a
                               % distância entre o centro da tela e esse
                               % ponto
    deg   = rad.*(180/pi); % Ângulo é *180/pi -- daria pra usar rad2deg
end
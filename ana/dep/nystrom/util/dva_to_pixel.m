function d = dva_to_pixel(data, options)
% Coverte dados em dva para pixel
% INPUT
% data:            dados de eye-tracking em pixels;
% dist:            distância entre o monitor e o sujeito, em cm. Padrão: 45;
% pixel_dva_ratio: quantos pixels são compreendidos por 1° (lembre: a relação
%                   não é linear, então para 2° não será o dobro);
% width:           largura da tela, em cm;
% res:             resolução da tela, em pixels (e.g. [1920, 1080] ou
%                   [1680, 1050], mas só importa a largura).
%
% OUTPUT:
%   d     - data converted from dva to pixels
%
% EXEMPLO: d2 = dva_to_pixel(data, 'dist', 45, 'width', 47.376, 'res', 1680);
%          d3 = dva_to_pixel(data, 'dist', 45, 'pixel_dva_ratio', 27.854);

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
        cpr = l/pdr;                       % cpr = cm_pixel_ratio (i.e., largura de 1 pixel em cm)
    
        d = dva2rad2pixel(data, cpr, options.dist);
    
    elseif ~isnan(options.width) && ~isnan(options.res(1))
    
        cpr = options.width/options.res(1); % cpr = cm_pixel_ratio 
        d = dva2rad2pixel(data, cpr, options.dist);
    
    else
        disp("Passe algum argumento adicional para transformar seus dados.")
        d = data;
    end
end

% Converte deg em radiano e depois em pixel
function pix = dva2rad2pixel(deg, c, dist)
    rad = tan(deg.*(pi/180));
    pix   = rad.*(dist/c);
end
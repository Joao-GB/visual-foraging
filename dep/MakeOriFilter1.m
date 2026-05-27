function [oriFilter, filtSize] = MakeOriFilter1(stimDiam_pix, aSigma, rSigma2, preferredOri)

    if nargin < 3,      rSigma2 = .76; end
    if nargin < 4, preferredOri = 90;  end

    filtSize = 2.^ceil(log2([stimDiam_pix, stimDiam_pix]));
    ImSize_x = filtSize(2); ImSize_y = filtSize(1);

    oriFilter = zeros(ImSize_y,ImSize_x);
    for id_X = 1:ImSize_x
        for id_Y = round(ImSize_y/2):ImSize_y

            x = id_X - ImSize_x/2; y = id_Y - ImSize_y/2;
            % Componente angular
            alpha = atan2d(y,x);
            angTerm = exp(-(alpha - preferredOri)^2 / aSigma^2);

            % COmponente radial
            if isnan(rSigma2)
                radTerm = 1;
            else
                r2 = x^2 + y^2; r2 = r2/(ImSize_x*ImSize_y);
                radTerm = exp(-r2 / rSigma2);
            end 
            fVal = angTerm * radTerm;

            oriFilter(id_Y,id_X) = fVal;
            oriFilter(ImSize_y-id_Y+1,ImSize_x-id_X+1) = fVal;
        end
    end
    oriFilter = rot90(oriFilter);
end
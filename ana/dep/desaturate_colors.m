function dim_colors = desaturate_colors(rgb_colors, saturation_factor, lightness_factor)
    % Desaturate RGB colors by reducing saturation in HSL space.
    % Automatically detects if input is normalized [0,1] or not [0,255].
    %
    % Inputs:
    %   rgb_colors: Nx3 matrix of RGB values (either 0-1 or 0-255)
    %   saturation_factor: Desaturation strength (0 = gray, 1 = no change)
    %
    % Output:
    %   dim_colors: Desaturated RGB in the same range as input.

    % Check if normalization is needed (if any value > 1)
    if any(rgb_colors(:) > 1)
        % Input is in 0-255 → normalize to 0-1
        rgb_normalized = double(rgb_colors) / 255;
        was_normalized = false;
    else
        % Input is already in 0-1
        rgb_normalized = rgb_colors;
        was_normalized = true;
    end

    % Initialize output
    dim_colors = zeros(size(rgb_colors));

    for i = 1:size(rgb_normalized, 1)
        R = rgb_normalized(i, 1);
        G = rgb_normalized(i, 2);
        B = rgb_normalized(i, 3);

        % Convert RGB to HSL
        [H, S, L] = rgb2hsl(R, G, B);

        % Apply saturation and lightness adjustments
        S = S * saturation_factor;
        if lightness_factor > 0
            L = L + (1 - L) * lightness_factor; % Move toward white
        elseif lightness_factor < 0
            L = L * (1 + lightness_factor);     % Move toward black
        end

        % Convert HSL back to RGB
        [R_new, G_new, B_new] = hsl2rgb(H, S, L);

        % Store result
        dim_colors(i, :) = [R_new, G_new, B_new];
    end

    % Revert to 0-255 if input was not normalized
    if ~was_normalized
        dim_colors = round(dim_colors * 255);
    end
end

function [H, S, L] = rgb2hsl(R, G, B)
    % Convert RGB [0-1] to HSL
    max_val = max([R, G, B]);
    min_val = min([R, G, B]);
    L = (max_val + min_val) / 2;

    if max_val == min_val
        H = 0;
        S = 0;
    else
        delta = max_val - min_val;
        if L > 0.5
            S = delta / (2 - max_val - min_val);
        else
            S = delta / (max_val + min_val);
        end

        if R == max_val
            H = (G - B) / delta + (G < B) * 6;
        elseif G == max_val
            H = (B - R) / delta + 2;
        else
            H = (R - G) / delta + 4;
        end
        H = H / 6;
    end
end

function [R, G, B] = hsl2rgb(H, S, L)
    % Convert HSL [0-1] to RGB
    if S == 0
        R = L;
        G = L;
        B = L;
    else
        if L < 0.5
            q = L * (1 + S);
        else
            q = L + S - L * S;
        end
        p = 2 * L - q;

        R = hue2rgb(p, q, H + 1/3);
        G = hue2rgb(p, q, H);
        B = hue2rgb(p, q, H - 1/3);
    end
end

function channel = hue2rgb(p, q, t)
    % Helper for hsl2rgb
    if t < 0
        t = t + 1;
    elseif t > 1
        t = t - 1;
    end

    if t < 1/6
        channel = p + (q - p) * 6 * t;
    elseif t < 1/2
        channel = q;
    elseif t < 2/3
        channel = p + (q - p) * (2/3 - t) * 6;
    else
        channel = p;
    end
end
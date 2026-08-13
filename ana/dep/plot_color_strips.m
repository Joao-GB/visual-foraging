function plot_color_strips(colors, strip_height)
    % Plot horizontal color strips for visualization
    % Inputs:
    %   colors: Nx3 matrix of RGB values (0-255 or 0-1)
    %   strip_height: Height of each color strip (default: 50 pixels)

    if nargin < 2
        strip_height = 50; % Default strip height
    end

    % Check if colors need normalization (if any value > 1)
    if any(colors(:) > 1)
        colors = double(colors) / 255; % Normalize to [0,1]
    end

    N = size(colors, 1);
    img = zeros(strip_height, N, 3); % Initialize image (height x width x RGB)

    % Assign each color to a vertical strip (repeated for height)
    for i = 1:N
        % Repeat the color strip_height times along the vertical axis
        img(:, i, 1) = colors(i, 1); % R channel
        img(:, i, 2) = colors(i, 2); % G channel
        img(:, i, 3) = colors(i, 3); % B channel
    end

    % Display the image
    figure;
    imshow(img);
    title('Color Strips Visualization');
    axis on;
    xticks(1:N);
    xlabel('Color Index');
    yticks([]);
end
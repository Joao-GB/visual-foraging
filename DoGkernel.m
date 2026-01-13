function kernel = DoGkernel(kernelSize, sigma_low, sigma_high)
    % kernelSize: size of kernel [height, width] or scalar for square
    % sigma_low: controls high-frequency cutoff (larger = lower frequencies pass)
    % sigma_high: controls low-frequency cutoff (smaller = higher frequencies pass)
    
    if numel(kernelSize) == 1
        kernelSize = [kernelSize, kernelSize];
    end
    
    % Create coordinate grid
    [x, y] = meshgrid(-(kernelSize(2)-1)/2:(kernelSize(2)-1)/2, ...
                      -(kernelSize(1)-1)/2:(kernelSize(1)-1)/2);
    
    % Create low-pass (wide Gaussian) and high-pass (narrow Gaussian)
    gauss_low = exp(-(x.^2 + y.^2) / (2 * sigma_low^2));
    gauss_high = exp(-(x.^2 + y.^2) / (2 * sigma_high^2));
    
    % Normalize
    gauss_low = gauss_low / sum(gauss_low(:));
    gauss_high = gauss_high / sum(gauss_high(:));
    
    % Band-pass = low-pass - high-pass
    kernel = gauss_low - gauss_high;
    
    % Normalize to zero mean (optional, but often desired for band-pass)
    kernel = kernel - mean(kernel(:));
    kernel = kernel / sum(abs(kernel(:))); % Normalize
end
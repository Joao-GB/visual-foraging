function filteredMat = butterFilter(mat, f_lo, f_hi, order)
    if nargin < 4, order = 4; end
    
    F = fft2(mat);
    Fshift = fftshift(F);
    
    [rows, cols] = size(mat);
    [u, v] = meshgrid( linspace(-0.5,0.5,cols), linspace(-0.5,0.5,rows) );
    
    freq = sqrt(u.^2 + v.^2);
    
    H_low  = 1 ./ (1 + (freq./f_hi).^(2*order));   % low-pass
    H_high = 1 ./ (1 + (f_lo./freq).^(2*order));   % high-pass
    H = H_low .* H_high;                           % band-pass
    
    H(freq==0) = 0;
    
    Ffiltered = Fshift .* H;
    
    filteredMat = real(ifft2(ifftshift(Ffiltered)));
end
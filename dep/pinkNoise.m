function pinkNoise = pinkNoise(N,M)
% Produz matriz NxM de ruído rosa
    % Gera ruído gaussiano
    rng('shuffle');
    whiteNoise = randn(N, M);
    
    % Gera a transformada de Fourier do sinal
    Nfft = fftshift(fft2(whiteNoise));
    
    % Constrói grid de frequência radial
    [fx, fy] = meshgrid( (-M/2:M/2-1)/M, (-N/2:N/2-1)/N );
    r = sqrt(fx.^2 + fy.^2);
    r(round(N/2)+1, round(M/2)+1) = 1;
    
    % Filtro de amplitude a ser aplicado na transformada
    Hpink = 1 ./ r;
    
    Npink_fft = Nfft .* Hpink;
    
    % Volta para o espaço original do sinal (por isso operações inversas)
    pinkNoise = real(ifft2(ifftshift(Npink_fft)));
end
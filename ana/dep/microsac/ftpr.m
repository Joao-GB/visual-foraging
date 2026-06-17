function xs = ftpr(x);
%
%
%
N = length(x);
if mod(N,2)~=1
    error('N is even in ftpr.m');
end
z = fftshift(fft(x));
phi = pi*(2*rand(1,N)-1);
z1 = z.*exp(i*phi);
z2re = real(z1+fliplr(z1))/2;
z2im = imag(z1-fliplr(z1))/2;
z2 = z2re + i*z2im;
xs = ifft(ifftshift(z2));

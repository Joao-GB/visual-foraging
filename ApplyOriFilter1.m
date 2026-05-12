function noiseIm_filt = ApplyOriFilter1(oriFilter, filtSize, noiseIm)

    meanSub = mean(noiseIm(:));
    noiseIm_fft = fftshift(fft2(noiseIm-meanSub,filtSize(1),filtSize(2)));    % fft and shift
    noiseIm_fft_filt = oriFilter .* noiseIm_fft;                                   % apply filter
    noiseIm_filt = real(ifft2(ifftshift(noiseIm_fft_filt)));                    % shift back
    noiseIm_filt = noiseIm_filt(1:size(noiseIm,1),1:size(noiseIm,2));
    st = std(noiseIm_filt(:));
%     sprintf('Pink std: %.4f', st)
    noiseIm_filt = (noiseIm_filt - mean(noiseIm_filt(:)))/st; 
end
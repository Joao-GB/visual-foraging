function noiseIm_filt = ApplyOriFilter(oriFilter, filtSize, noiseIm)

    meanSub = mean(noiseIm(:));
    noiseIm_fft = fftshift(fft2(noiseIm-meanSub,filtSize(1),filtSize(2)));    % fft and shift
    noiseIm_fft_filt = oriFilter .* noiseIm_fft;                                   % apply filter
    noiseIm_filt = real(ifft2(ifftshift(noiseIm_fft_filt)));                    % shift back
    noiseIm_filt = noiseIm_filt(1:size(noiseIm,1),1:size(noiseIm,2));
    noiseIm_filt = noiseIm_filt+meanSub;
    noiseIm_filt = noiseIm_filt-min(noiseIm_filt(:));
    noiseIm_filt = noiseIm_filt./max(noiseIm_filt(:));
    
end
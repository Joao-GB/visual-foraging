sec_dev_den_x = diff(diff(x_denoised));
acc_filt_x = diff(diff(st.filt_data(1,:)));

T = 1:17;
L = length(T);
filt_counts = zeros(1,L);
tvd_counts = zeros(1,L);
for i=1:L
    [~, tvd_peaks] = findpeaks(sec_dev_den_x(abs(sec_dev_den_x) > 10^(-T(i))));

    tvd_counts(i)  = length(tvd_peaks);
    [~, filt_peaks] = findpeaks(acc_filt_x(abs(acc_filt_x) > 10^(-T(i))));
    filt_counts(i) = length(filt_peaks);

end
figure;
plot(T, filt_counts);
hold on;
plot(T, tvd_counts);
legend('Filt data acc','TV denoising acc')
title('Number of acc peaks above 10^{(-x)}')
xlabel('x (exponent)')
ylabel('Counts')
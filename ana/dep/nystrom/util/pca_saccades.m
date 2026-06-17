function score = pca_saccades(data, fs)

    w = 50;
    order = 2;
    filt_size = .01*fs + 1;

    x = double(data(1,:));
    y = double(data(2,:));
    r = sqrt(x.^2 + y.^2);

    [~,~,dx] = unsharp_masking(x);
    [~,~,dy] = unsharp_masking(y);

    mmx = movmean(x, w);
    mmy = movmean(y, w);
    mmr = movmean(r, w);

    mvx = movvar(x, w);
    mvy = movvar(y, w);
    mvr = movvar(r, w);

%     dat(1, :) = filtfilt(b, 1, x);
%     dat(2, :) = filtfilt(b, 1, y);

%     vx = filtfilt(gSG1, 1, x);
%     vy = filtfilt(gSG1, 1, y);
% 
%     ax = filtfilt(gSG2, 1, x);
%     ay = filtfilt(gSG2, 1, y);

    [v, a] = nystrom_vel_n_acc(data, fs,order, filt_size, 'all vel and acc, please!');
    vel = v. vel; vx = v.vx; vy = v.vy;
    acc = a.acc; ax = a.ax; ay = a.ay;
%     features = [x;y;r;vx;vy;vel;ax;ay;acc].';
%     features = [vx;vy;vel;ax;ay;acc].';
    features = [mmx;mmy;mmr;mvx;mvy;mvr;vx;vy;vel;ax;ay;acc].';

    z_feat = zscore(features);

    [~, score, ~] = pca(z_feat);

    % Plot first 2-3 PCs
    figure;
    scatter(score(:,1), score(:,2), 5, 'filled');
    xlabel('PC1'); ylabel('PC2');
    title('PCA of eye-tracking 14 features');
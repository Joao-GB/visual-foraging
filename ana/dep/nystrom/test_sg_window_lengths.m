function test_sg_window_lengths(data, fs, order, len_values)

    figure; hold on;
    stats = [];

    for len = len_values
        if mod(len, 2) == 0
            warning("Skipping even SG frame length: %d", len);
            continue
        end

        % Compute velocity using your SG method
        [vel, ~] = nystrom_vel_n_acc(data, fs, order, len);

        % Run adaptive thresholding
        v_peak = nystrom_v_peak(vel, 120, 1);  % or another init threshold

        % Get fixation-only portion for stats
        v_below = vel(vel < v_peak);
        med_val = median(v_below);
        mad_val = 1.4826 * mad(v_below);

        % Save stats
        stats = [stats; len, v_peak, med_val, mad_val];

        % Plot histogram
        subplot(ceil(length(len_values)/2), 2, find(len_values==len)); 
        histogram(vel, 100);
        xline(v_peak, 'r', 'LineWidth', 1.5);
        title(sprintf('Len=%d | v_peak=%.1f | MAD=%.1f', len, v_peak, mad_val));
        xlabel('Velocity (°/s)'); ylabel('Count');
    end

    % Summary plot
    figure;
    plot(stats(:,1), stats(:,2), '-o', 'LineWidth', 2); hold on;
    plot(stats(:,1), stats(:,4), '-x', 'LineWidth', 2);
    legend('Final Threshold (v\_peak)', 'MAD (below threshold)');
    xlabel('SG Frame Length'); ylabel('Value (°/s)');
    title('Effect of SG Length on Adaptive Thresholding');
    grid on;

end

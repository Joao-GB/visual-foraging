function [v_peaks, amp] = get_v_peak_n_amp(data, vel, sac_lims, mode)
if nargin < 4
    mode = 'peaks';
end
    L = size(sac_lims,2);
    amp = zeros(1, L);
    v_peak = zeros(1, L);
    vp_idx = zeros(1, L);
    for i = 1:L
        amp(i) = sqrt( (data(1,sac_lims(2,i))-data(1,sac_lims(1,i)))^2 + (data(2,sac_lims(2,i))-data(2,sac_lims(1,i)))^2 );
        if strcmp(mode, 'peaks')
            [v_peak(i), vp_idx(i)] = findpeaks(vel(sac_lims(1,i):sac_lims(2,i)), SortStr='descend', NPeaks=1);
        else
            [v_peak(i), vp_idx(i)] = max((vel(sac_lims(1,i):sac_lims(2,i))));
        end

        vp_idx(i) = vp_idx(i) + sac_lims(1,i) - 1;
    end
    v_peaks.v_peak = v_peak;
    v_peaks.idx    = vp_idx;
end
    
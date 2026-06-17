function [fix, non_fix] = nystrom_get_fixations(data, nan_idx, min_fix)
if nargin <=2
    min_fix = 1;
end
[nan_on, nan_off] = get_start_stop(nan_idx);
nan_lims = [nan_on; nan_off];

all_itv = [1;length(data(1,:))];
f_lims = trim_overlapping_events(all_itv, nan_lims, 1);

f_on = f_lims(1,:);
f_off = f_lims(2,:);
L = length(data);

non_fix = [0;0];
for i = 1:length(f_on)
    if f_off(i) - f_on(i) < min_fix
        if f_on(i) ~= 1;  f_on(i) = f_on(i) + 1;   end
        if f_off(i) ~= L; f_off(i) = f_off(i) - 1; end
        non_fix = [non_fix(1,:) f_on(i); non_fix(2,:) f_off(i)];
        f_on(i)  = NaN;
        f_off(i) = NaN;
    end
end
non_fix(:,1) = [];
f_on(isnan(f_on))   = [];
f_off(isnan(f_off)) = [];

fix = [f_on; f_off];
function pur = get_pursuits(data, fix, thr)
    pur = [];
    for i=1:length(fix(1,:))
        [diam, ~] = trajectory_dispersion(data(:, fix(1,i):fix(2,i)));
        if diam > thr
            pur = [pur fix(:, i)];
        end
    end


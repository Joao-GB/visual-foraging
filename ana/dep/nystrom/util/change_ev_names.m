function ev = change_ev_names(ev)
    ev_possible_names = {{'fixations', 'fix', 'fixs', 'fixation', 'fixational movs', 'fixational mov', 'fixational movement', 'fixational movements'},...
                         {'saccades', 'sac', 'sacs', 'saccade', 'saccadic', 'saccadic mov', 'saccadic movs', 'saccadic movement', 'saccadic movements'},...
                         {'blinks', 'blk', 'blks', 'blink', 'noise', 'blink and noise', 'blinks and noise'},...
                         {'pursuits', 'pur', 'purs', 'pursuit', 'smooth pursuit', 'smooth pursuits'},...
                         {'psos', 'pso', 'glissade', 'glissades', 'oscillation', 'oscillations', 'post-saccadic oscillation', 'post-saccadic oscillations', 'post saccadic oscillation', 'post saccadic oscillations'}};
    for i = 1:length(ev)
        for j = 1:length(ev_possible_names)
            if any(strcmp(ev_possible_names{j}, ev{i}))
                aux = ev_possible_names{j};
                ev{i} = aux{1};
                break;
            end
        end
    end
end
function nums = foragingSearchSubj(outFolder, params)

    prefix = [params.matPreffix ];
    files = dir(fullfile(outFolder, [prefix '*']));

    nums = [];

    for k = 1:numel(files)
        name = files(k).name;

        % Extract the number immediately after the prefix
        tok = regexp(name, ['^' regexptranslate('escape', prefix) '(\d+)'], ...
                     'tokens', 'once');

        if ~isempty(tok)
            nums(end+1) = str2double(tok{1}); %#ok<AGROW> 
        end
    end

    nums = sort(unique(nums));
end
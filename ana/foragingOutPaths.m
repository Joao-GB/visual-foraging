function [matPath, edfPath] = foragingOutPaths(folder, subject)

% Procura por arquivos cujo nome começa com
%   s## e retorna matPath
%   e## e retorna edfPath
%
% Behavior:
% - If none found: prints message and returns []
% - If one found: auto-selects
% - If multiple: lists and prompts user to choose
%
% Inputs:
%   folder  - path to directory (char or string)
%   subject - numeric or string (e.g., 0, 1, '03')
%
% Outputs:
%   m - full path to selected s## file (or [])
%   e - full path to selected e## file (or [])

    if isnumeric(subject)
        subject_str = sprintf('%02d', subject);
    else
        subject_str = sprintf('%02d', str2double(subject));
    end

    s_pattern = ['s' subject_str];
    e_pattern = ['e' subject_str];

    files = dir(folder);
    files = files(~[files.isdir]); % ignore folders

    names = {files.name};

    s_idx = startsWith(names, s_pattern);
    e_idx = startsWith(names, e_pattern);

    s_files = files(s_idx);
    e_files = files(e_idx);

    % --- helper function ---
    function selected = handleGroup(fileStruct, label)
        if isempty(fileStruct)
            fprintf('Nenhum arquivo encontrado para o sujeito %s (%s)\n', subject_str, label);
            selected = [];
            return;
        end

        if numel(fileStruct) == 1
            fprintf('%s: selecionado automaticamente -> %s\n', label, fileStruct(1).name);
            selected = fullfile(fileStruct(1).folder, fileStruct(1).name);
            return;
        end

        fprintf('\nMais de um arquivo encontrado para %s (sujeito %s):\n', label, subject_str);
        for i = 1:numel(fileStruct)
            fprintf('[%d] %s\n', i, fileStruct(i).name);
        end

        while true
            idx = input(sprintf('Escolha o índice do arquivo %s: ', label));
            if isnumeric(idx) && isscalar(idx) && idx >= 1 && idx <= numel(fileStruct)
                selected = fullfile(fileStruct(idx).folder, fileStruct(idx).name);
                break;
            else
                fprintf('Entrada inválida, tente novamente.\n');
            end
        end
    end

    % --- process both groups ---
    matPath = handleGroup(s_files, 'm (s-files)');
    edfPath = handleGroup(e_files, 'e (e-files)');
end
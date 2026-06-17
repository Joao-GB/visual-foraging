function filesSelected = foragingFindFiles(folder, st, ext)
% FORAGINGFINDFILES Find the highest common-version set of files.
%
% filesSelected = foragingFindFiles(folder, st, ext)
%
% Searches a folder for files matching each pair
% {st{k}, ext{k}}. The selected files must all share the
% same version number '_vX'.
%
% Priority order:
%   1) Files containing neither '_r' nor '_i'
%   2) Files containing '_r'
%   3) Files containing '_i'
%
% Within the highest-priority non-empty group, the function
% selects the largest version number that exists for ALL
% requested file types.
%
% Inputs
% -------
% folder : char/string
%     Folder to search.
%
% st : cell array of strings
%     Required substrings in filenames.
%
% ext : cell array of strings
%     Required file extensions corresponding to st.
%
% Output
% -------
% filesSelected : cell array
%     Full paths of the selected files.
%
% Example
% -------
% files = foragingFindFiles( ...
%     'C:\data', ...
%     {'data','config'}, ...
%     {'.mat','.txt'});

    if numel(st) ~= numel(ext)
        error('st and ext must have the same length.');
    end

    D = dir(folder);
    D = D(~[D.isdir]);

    names = {D.name};

    filesSelected = {};

    % Priority order
    for priority = 1:3

        versionSets = cell(1,numel(st));
        fileMaps = cell(1,numel(st));

        validPriority = true;

        for k = 1:numel(st)

            versionSets{k} = [];
            fileMaps{k} = containers.Map('KeyType','double', ...
                                         'ValueType','char');

            foundAnyMatch = false;

            for i = 1:numel(names)

                fname = names{i};

                if ~contains(fname, st{k})
                    continue
                end

                [~,~,fext] = fileparts(fname);

                if ~strcmpi(fext, ext{k})
                    continue
                end

                foundAnyMatch = true;

                hasR = contains(fname,'_r');
                hasI = contains(fname,'_i');

                switch priority
                    case 1
                        if hasR || hasI
                            continue
                        end

                    case 2
                        if ~hasR
                            continue
                        end

                    case 3
                        if ~hasI
                            continue
                        end
                end

                tok = regexp(fname,'_v(\d+)','tokens','once');

                if isempty(tok)
                    version = 0;
                else
                    version = str2double(tok{1});
                end
                
                versionSets{k}(end+1) = version;
                
                % Keep highest occurrence if duplicated
                fileMaps{k}(version) = fname;
            end

            if ~foundAnyMatch
                error('No file found containing "%s".', st{k});
            end
        end

        if ~validPriority
            continue
        end

        commonVersions = versionSets{1};

        for k = 2:numel(versionSets)
            commonVersions = intersect(commonVersions, ...
                                       versionSets{k});
        end

        if isempty(commonVersions)
            continue
        end

        bestVersion = max(commonVersions);

        filesSelected = cell(1,numel(st));

        for k = 1:numel(st)
            filesSelected{k} = ...
                fullfile(folder, fileMaps{k}(bestVersion));
        end

        return
    end

    warning('No common version found for all requested files.');
end
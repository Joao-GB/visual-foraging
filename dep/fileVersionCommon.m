function suffix = fileVersionCommon(baseNames, extensions)
    % baseNames:  cell array de strings com os nomes base (ex: {'sc_Sub1_cond1', 'sc_Sub1_All'})
    % extensions: cell array ou string única com as extensões (ex: '.png' ou {'.png', '.png'})
    
    % Se passarem apenas uma extensão em string, replica para todos
    if ischar(extensions) || isstring(extensions)
        extensions = repmat({char(extensions)}, size(baseNames));
    end
    
    numFiles = numel(baseNames);
    
    % Passo 1: Verificar se a versão "master" (sem sufixo) está livre para TODOS
    allFree = true;
    for i = 1:numFiles
        if exist([baseNames{i}, extensions{i}], 'file')
            allFree = false;
            break;
        end
    end
    
    % Se nenhum arquivo existe ainda, o sufixo é vazio (versão original)
    if allFree
        suffix = '';
        return;
    end
    
    % Passo 2: Se algum já existe, procura a primeira versão '_vX' comum livre
    counter = 1;
    while true
        suffixCandidate = sprintf('_v%d', counter);
        allFree = true;
        
        for i = 1:numFiles
            checkName = [baseNames{i}, suffixCandidate, extensions{i}];
            if exist(checkName, 'file')
                allFree = false;
                break; % Esse counter já está ocupado por alguém, pula para o próximo
            end
        end
        
        if allFree
            suffix = suffixCandidate;
            return;
        end
        
        counter = counter + 1;
    end
end
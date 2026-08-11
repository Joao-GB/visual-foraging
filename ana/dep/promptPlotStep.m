function [skipMode, closePlots] = promptPlotStep(sectionName)
    closePlots = false;
    
    % Main prompt
    choice = questdlg(sprintf('Exibir: %s?', sectionName), ...
        'Gerenciamento de Plots', ...
        'Prosseguir', 'Pular próximo', 'Pular todos', 'Prosseguir');
    
    switch choice
        case 'Prosseguir'
            skipMode = 'none';
            % Secondary prompt only if user proceeds
            closeChoice = questdlg('Fechar plots abertos atualmente?', ...
                'Gerenciamento de Plots', ...
                'Sim', 'Não', 'Não');
            if strcmp(closeChoice, 'Sim')
                closePlots = true;
            end
        case 'Pular próximo'
            skipMode = 'this';
        case 'Pular todos'
            skipMode = 'all';
        otherwise % Handle window 'X' close button
            skipMode = 'this';
    end
end
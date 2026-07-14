function renderPSAforagingNumSeen(barMatrix, ~, uniqueSeen, ylabelStr, figName, plotTitle, drP, isSens)
    figure('Name', figName, 'Color', 'w', 'Position', [100 200 800 500]);
    
    if isSens
        b = bar(barMatrix, 0.3, 'grouped', 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
    else
        b = bar(barMatrix, 'grouped', 'EdgeColor', [0 0 0]);
    end
    
    % Configuração de paleta institucional
    b(1).FaceColor = drP.darkBlue;   % Primeira barra do par: Sacádico
    b(2).FaceColor = drP.paleBrown;  % Segunda barra do par: Não-sacádico
    
    % Estética básica dos eixos
    set(gca, 'TickDir', 'out', 'Box', 'off');
    grid on;
    
    set(gca, 'XTick', 1:length(uniqueSeen));
    xticklabels(string(uniqueSeen));
    
    xlabel('Número de estímulos vistos');
    ylabel(ylabelStr);
    title(plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
    legend({'Sacádico', 'Não-sacádico'}, 'Location', 'northeast', 'Box', 'off');
    
    if isSens
        localMax = max(barMatrix(:));
        if localMax > 0
            ylim([0 localMax * 1.2]);
        else
            ylim([0 1]);
        end
    else
        ylim([0 100]);
    end
end
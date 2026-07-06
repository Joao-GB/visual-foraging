function PSA = plotPSAmain(trl, drP)
    titlePSA  = 'Efeito pré-sacádico em tarefa de forrageamento';
    xlabelPSA = ["Forrag.", "Sacádico", "Não-sacádico"];
    ylabelPSA = 'Acertos (%)';
    [main, counts] = getPSAeffect(trl);
    PSA.main = main; PSA.main.counts = counts;
    
    figure('Name', 'PSA main effect', 'Color', 'w', 'Position', [100 200 1000 450]);
    t = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'normal');
    
    % Subplot 1: acertos
    ax1 = nexttile;
    barValues = (counts(1,:)./counts(2,:)) * 100;
    
    b1 = bar(barValues, 'FaceColor', 'flat');
    b1.CData(1,:) = [1 1 1];
    b1.CData(2,:) = drP.darkBlue;
    b1.CData(3,:) = drP.paleBrown;
    b1.EdgeColor  = [0 0 0];
    
    set(ax1, 'TickDir', 'out', 'Box', 'off')
    xticklabels(xlabelPSA);
    ylabel(ylabelPSA); ylim([0 100]);
    title('Acurácia');
    grid on;

    % Subplot 2: d-prime
    ax2 = nexttile;
%     dPrimeValues = [main.for.d, main.sacc.d, main.nSacc.d];
    dPrimeValues = [main.sacc.d, main.nSacc.d];
    
    b2 = bar(dPrimeValues, .3, 'FaceColor', 'flat', 'LineWidth', 1.5);
%     b2.CData(1,:) = [1 1 1];
%     b2.CData(2,:) = drP.darkBlue;
%     b2.CData(3,:) = drP.paleBrown;
    b2.CData(1,:) = drP.darkBlue;
    b2.CData(2,:) = drP.paleBrown;
    b2.EdgeColor  = [0 0 0];
    
    set(ax2, 'TickDir', 'out', 'Box', 'off')
%     xticklabels(["Forrag.", "Sacádico", "Não-sacádico"]);
    xticklabels(["Sacádico", "Não-sacádico"]);
    ylabel('Sensibilidade (d'')'); 
    
    if max(dPrimeValues) > 0
        ylim([0 max(dPrimeValues) * 1.2]);
    end
    title('Sensibilidade Perceptual');
    grid on;

    title(t, titlePSA, 'FontSize', 14, 'FontWeight', 'bold');
end
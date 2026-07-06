function plotPSAforagingPerformance(trl, PSA, drP)
% esperado que, se errou forrageamento, talvez não estivesse prestando
% atenção então também erre os demais?

    correctFor = logical(PSA.main.for.idx{1,1} + PSA.main.for.idx{2,2});
    incorrectFor = ~correctFor;
    cFtrl = trl(correctFor);
    iFtrl = trl(incorrectFor);
    [~, cFcounts] = getPSAeffect(cFtrl);
    [~, iFcounts] = getPSAeffect(iFtrl);
    
    titlePSA  = 'Efeito pré-sacádico e desempenho no forrageamento';
    xlabelPSA = ["Sacádico", "Não-sacádico"];
    ylabelPSA = 'Acertos (%)';
    
    figure('Name', 'PSA foraging effect', 'Color', 'w', 'Position', [100 200 900 450]);
    
    t = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'normal');
    
    % Subplot 1: forrageamento completo
    ax1 = nexttile;
    % Apenas as colunas correspondentes de cFcounts
    barValues_cF = (cFcounts(1,2:3)./cFcounts(2,2:3)) * 100;
    
    b1 = bar(barValues_cF, 'FaceColor', 'flat');
    b1.CData(1,:) = drP.darkBlue;
    b1.CData(2,:) = drP.paleBrown;
    
    set(ax1, 'TickDir', 'out', 'Box', 'off')
    xticklabels(xlabelPSA);
    ylabel(ylabelPSA); ylim([0 100]);
    title('Forrageamento Correto');
    grid on;

    % Subplot 2: forrageamento incompleto
    ax2 = nexttile;
    barValues_iF = (iFcounts(1,2:3)./iFcounts(2,2:3)) * 100;
    
    b2 = bar(barValues_iF, 'FaceColor', 'flat');
    b2.CData(1,:) = drP.darkBlue;
    b2.CData(2,:) = drP.paleBrown;
    
    set(ax2, 'TickDir', 'out', 'Box', 'off')
    xticklabels(xlabelPSA);
    ylim([0 100]);
    title('Forrageamento Incorreto');
    grid on;
    
    % Título global
    title(t, titlePSA, 'FontSize', 13, 'FontWeight', 'bold');
    
end
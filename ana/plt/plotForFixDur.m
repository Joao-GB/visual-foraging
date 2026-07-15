function plotForFixDur(trl, mat)
    histLen    = {trl.forHistLen};
    histFixDur = {trl.forHistFixDur};
    histCat    = {trl.forHistCat};
    histFixDurNCat = cellfun(@(x,y,z)[x(1:z); y(1:z)],histFixDur, histCat, histLen, 'UniformOutput',false);
    histFixDurNCat = [histFixDurNCat{:}];

    tgtDurDistr   = histFixDurNCat(1,histFixDurNCat(2,:) == 1);
    noTgtDurDistr = histFixDurNCat(1,histFixDurNCat(2,:) == 0);

    forFixDurNAns = [[trl.forProbeFixDur]; [trl.forProbeHit]];

    hitDurDistr  = forFixDurNAns(1, forFixDurNAns(2,:) == 1);
    missDurDistr = forFixDurNAns(1, forFixDurNAns(2,:) == 0);
    

    if isfield(mat, 'drP')
        cTarget = mat.drP.darkBlue;
        cDist   = mat.drP.paleBrown;
    else
        cTarget = [0.00, 0.45, 0.74];
        cDist   = [0.80, 0.53, 0.22];
    end
    
    figure('Name', 'Foraging fixation distribution effects', 'Color', 'w', 'Position', [150, 150, 950, 450]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');

    nexttile;
    drawForagingBoxplot(tgtDurDistr, noTgtDurDistr, {'Alvo', 'Distrator'}, ...
        'Duração por categoria do estímulo', cTarget, cDist);

    nexttile;
    drawForagingBoxplot(hitDurDistr, missDurDistr, {'Acerto', 'Erro'}, ...
        'Duração por desempenho', cTarget, cDist);
end

function drawForagingBoxplot(distr1, distr2, boxLabels, plotTitle, col1, col2)
    data = [distr1(:); distr2(:)];
    groups = [ones(numel(distr1), 1); 2 * ones(numel(distr2), 1)];
    
    h = boxplot(data, groups, ...
        'Colors', [col1; col2], ...
        'Labels', boxLabels, ...
        'Symbol', 'o', ...
        'OutlierSize', 4, ...
        'Widths', 0.5);
    
    set(h, 'LineWidth', 1.5);
    set(gca, 'TickDir', 'out', 'Box', 'off');
    ylabel('Duração da Fixação');
    title(plotTitle, 'FontSize', 11, 'FontWeight', 'bold');
    
    yMax = max(data);
    yRange = max(data) - min(data);
    if isempty(yRange) || yRange == 0, yRange = 1; end
    ylim([0, yMax + (yRange * 0.15)]);
    
    text(1, yMax + (yRange * 0.05), sprintf('n=%d', numel(distr1)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', 'Color', col1);
    
    text(2, yMax + (yRange * 0.05), sprintf('n=%d', numel(distr2)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', 'Color', col2);
end
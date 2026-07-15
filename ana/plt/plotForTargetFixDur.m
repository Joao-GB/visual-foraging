function plotForTargetFixDur(trl, mat)
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
    
    data = [tgtDurDistr(:); noTgtDurDistr(:)];
    
    % Vetor de agrupamento: 1 para tgt, 2 para distrator
    groups = [ones(numel(tgtDurDistr), 1); 2 * ones(numel(noTgtDurDistr), 1)];
    
    figure('Name', 'Foraging category-fixation effect', 'Color', 'w', 'Position', [200, 200, 500, 450]);
    ax = axes();
    h = boxplot(data, groups, ...
        'Colors', [cTarget; cDist], ...
        'Labels', {'Alvo', 'Distrator'}, ...
        'Symbol', 'o', ...
        'OutlierSize', 4, ...
        'Widths', 0.5);
    
    set(h, 'LineWidth', 1.5);
    
    set(ax, 'TickDir', 'out', 'Box', 'off');
    ylabel('Duração da Fixação (s)');
    title('Duração da fixação no forrageamento', 'FontSize', 12, 'FontWeight', 'bold');
    
    yMax = max(data);
    yRange = max(data) - min(data);
    ylim([0, yMax + (yRange * 0.15)]); % Dá 15% de respiro no topo
    
    text(1, yMax + (yRange * 0.05), sprintf('n=%d', numel(tgtDurDistr)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', 'Color', cTarget);
        
    text(2, yMax + (yRange * 0.05), sprintf('n=%d', numel(noTgtDurDistr)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', 'Color', cDist);
end
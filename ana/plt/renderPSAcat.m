function renderPSAcat(plotData, ylabelStr, figName, plotTitle, drP, isSens)
    figure('Name', figName, 'Color', 'w', 'Position', [50 50 700 700]);
    
    ax_core_right = []; ax_marg_right = []; ax_core_bottom = []; ax_marg_bottom = [];

    allVals = cellfun(@(x) max(x(:)), plotData, 'UniformOutput',false);
    allVals = [allVals{:}];
    globalMax = max(allVals(:));
    if isnan(globalMax) || globalMax == 0, globalMax = 1; end

    t = tiledlayout(10, 10, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    rowLabels = {'P: alvo', 'P: distrator', 'Marginal (S)'};
    colLabels = {'S: alvo', 'S: distrator', 'Marginal (P)'};
    
    rowStarts = [1, 4, 8];
    colStarts = [1, 4, 8];
    
    for r = 1:3
        for c = 1:3
            if r == 3 && c == 3, continue; end
            
            tileIdx = (rowStarts(r) - 1) * 10 + colStarts(c);
            ax = nexttile(tileIdx, [3, 3]); 
            
            if r == 1 && c == 2, ax_core_right = ax;  end
            if r == 1 && c == 3, ax_marg_right = ax;  end
            if r == 2 && c == 1, ax_core_bottom = ax; end
            if r == 3 && c == 1, ax_marg_bottom = ax; end
            
            % Ajusta espessura baseado no tipo de plot
            if isSens
                b = bar(plotData{r,c}, 0.3, 'FaceColor', 'flat', 'LineWidth', 1.5);
            else
                b = bar(plotData{r,c}, 'FaceColor', 'flat');
            end
            
            b.EdgeColor = [0 0 0];
            b.CData(1,:) = drP.darkBlue;
            b.CData(2,:) = drP.paleBrown;
            set(ax, 'TickDir', 'out', 'Box', 'off');
            xticklabels({'S', 'N'});
            
            % Escalonamento dinâmico do eixo Y
            if isSens
                ylim([0 globalMax * 1.2]);
            else
                ylim([0 115]);
            end
            
            grid on;
            if c == 1, ylabel({['\bf{', rowLabels{r}, '}\rm'], ylabelStr}, 'Interpreter', 'tex'); end
            if r == 1, title(colLabels{c}); end
            if r == 3 || c == 3, set(ax, 'Color', [0.93 0.93 0.95]); end
        end
    end
    
    title(t, plotTitle, 'FontSize', 14, 'FontWeight', 'bold');
    
    % Linhas para separar as marginais
    drawnow; 
    tPos = t.Position;
    
    % Linha Vertical 
    lineX = ax_core_right.Position(1) + ax_core_right.Position(3) + (ax_marg_right.Position(1) - (ax_core_right.Position(1) + ax_core_right.Position(3)))/2;
    annotation('line', [lineX lineX], [tPos(2) tPos(2)+tPos(4)], ...
        'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.2);
        
    % Linha Horizontal
    row2_bottom = ax_core_bottom.Position(2);
    row3_top    = ax_marg_bottom.Position(2) + ax_marg_bottom.Position(4);
    lineY = row3_top + (row2_bottom - row3_top)/2;
    
    annotation('line', [tPos(1) tPos(1)+tPos(3)], [lineY lineY], ...
        'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.2);
end
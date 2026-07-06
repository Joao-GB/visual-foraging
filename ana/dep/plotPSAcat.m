function plotPSAcat(fixCat, sCat, sHit, nHit, drP)
    fixCat = fixCat(:); sCat = sCat(:); sHit = sHit(:); nHit = nHit(:);
    
    plotData = cell(3,3);
    
    % Calcula o plot principal 2x2
    for r = 1:2 % linhas: categoria do pré-probe (1 = alvo, 2 = distrator)
        fVal = 2 - r; % mapeia 1 -> fixCat=1, 2 -> fixCat=0
        
        for c = 1:2 % colunas: categoria do probe (1 = alvo, 2 = distrator)
            sVal = 2 - c; % 1 -> sCat=1, 2 -> sCat=0
            
            mask = (fixCat == fVal & sCat == sVal);
            
            if any(mask)
                % Calcula porcentagem de hits
                pct_s = (sum(sHit(mask)) / sum(mask)) * 100;
                pct_n = (sum(nHit(mask)) / sum(mask)) * 100;
                plotData{r,c} = [pct_s, pct_n];
            else
                plotData{r,c} = [0, 0];
            end
        end
    end
    
    % Calcula marginais das colunas (i.e., a última linha)
    for c = 1:2
        sVal = 2 - c;
        mask = (sCat == sVal);
        if any(mask)
            plotData{3,c} = [(sum(sHit(mask))/sum(mask))*100, (sum(nHit(mask))/sum(mask))*100];
        else
            plotData{3,c} = [0, 0];
        end
    end
    
    % Calcula marginais das linhas (i.e., última coluna)
    for r = 1:2
        fVal = 2 - r;
        mask = (fixCat == fVal);
        if any(mask)
            plotData{r,3} = [(sum(sHit(mask))/sum(mask))*100, (sum(nHit(mask))/sum(mask))*100];
        else
            plotData{r,3} = [0, 0];
        end
    end

    % Plota o 3x3 com o esquema das 10 colunas
    figure('Name', 'PSA category effect', 'Color', 'w', 'Position', [50 50 700 700]);
    
    % Create variables to store references
    ax_core_right = []; ax_marg_right = []; ax_core_bottom = []; ax_marg_bottom = [];

    % 1. CORRECTION: Define a 10x10 structural grid
    t = tiledlayout(10, 10, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    rowLabels = {'P: alvo', 'P: distrator', 'Marginal (S)'};
    colLabels = {'S: alvo', 'S: distrator', 'Marginal (P)'};
    
    % 2. CORRECTION: Define both row and column starts to create a symmetric grid
    rowStarts = [1, 4, 8];
    colStarts = [1, 4, 8];
    
    for r = 1:3
        for c = 1:3
            if r == 3 && c == 3, continue; end
            
            tileIdx = (rowStarts(r) - 1) * 10 + colStarts(c);
            ax = nexttile(tileIdx, [3, 3]); 
            
            % --- Keep the rest of your inner loop code exactly the same ---
            if r == 1 && c == 2, ax_core_right = ax;  end
            if r == 1 && c == 3, ax_marg_right = ax;  end
            if r == 2 && c == 1, ax_core_bottom = ax; end
            if r == 3 && c == 1, ax_marg_bottom = ax; end
            
            b = bar(plotData{r,c}, 'FaceColor', 'flat');
            b.EdgeColor = [0 0 0];
            b.CData(1,:) = drP.darkBlue;
            b.CData(2,:) = drP.paleBrown;
            set(ax, 'TickDir', 'out', 'Box', 'off');
            xticklabels({'S', 'N'});
            ylim([0 100]);
            grid on;
            if c == 1, ylabel({['\bf{', rowLabels{r}, '}\rm'], 'Acertos (%)'}, 'Interpreter', 'tex'); end
            if r == 1, title(colLabels{c}); end
            if r == 3 || c == 3, set(ax, 'Color', [0.93 0.93 0.95]); end
        end
    end
    
    title(t, 'Efeito de categoria dos estímulos', 'FontSize', 14, 'FontWeight', 'bold');
    
    %% --- 5. DRAW SEPARATING LINES FOR MARGINALS ---
    drawnow; 
    
    tPos = t.Position;
    
    % A. Vertical Line (Dead-center between Col 2 right-edge and Col 3 left-edge)
    lineX = ax_core_right.Position(1) + ax_core_right.Position(3) + (ax_marg_right.Position(1) - (ax_core_right.Position(1) + ax_core_right.Position(3)))/2;
    annotation('line', [lineX lineX], [tPos(2) tPos(2)+tPos(4)], ...
        'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.2);
        
    % B. CORRECTION: Horizontal Line (Dead-center between Row 2 bottom-edge and Row 3 top-edge)
    % Note: MATLAB's Y-axis positions start from the BOTTOM of the figure canvas.
    row2_bottom = ax_core_bottom.Position(2);
    row3_top    = ax_marg_bottom.Position(2) + ax_marg_bottom.Position(4);
    lineY = row3_top + (row2_bottom - row3_top)/2;
    
    annotation('line', [tPos(1) tPos(1)+tPos(3)], [lineY lineY], ...
        'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.2);
end
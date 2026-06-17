function getPSAcat(fixCat, sCat, sHit, nHit)
    % Ensure inputs are column vectors
    fixCat = fixCat(:); sCat = sCat(:); sHit = sHit(:); nHit = nHit(:);
    
    % Initialize a 3x3 cell array to hold data for each plot: [sHit_pct, nHit_pct]
    plotData = cell(3,3);
    
    % --- 1. CALCULATE THE CORE 2x2 CONDITIONS ---
    for r = 1:2 % Rows: Fixation Category (1 = Target, 2 = Distractor)
        fVal = 2 - r; % Maps row 1 -> fixCat=1, row 2 -> fixCat=0
        
        for c = 1:2 % Columns: Saccade Category (1 = Target, 2 = Distractor)
            sVal = 2 - c; % Maps col 1 -> sCat=1, col 2 -> sCat=0
            
            % Logical mask for current combination
            mask = (fixCat == fVal & sCat == sVal);
            
            if any(mask)
                % Calculate percentage of hits for both types
                pct_s = (sum(sHit(mask)) / sum(mask)) * 100;
                pct_n = (sum(nHit(mask)) / sum(mask)) * 100;
                plotData{r,c} = [pct_s, pct_n];
            else
                plotData{r,c} = [0, 0];
            end
        end
    end
    
    % --- 2. CALCULATE COLUMN MARGINALS (Bottom Row) ---
    for c = 1:2
        sVal = 2 - c;
        mask = (sCat == sVal); % Ignores fixation category
        if any(mask)
            plotData{3,c} = [(sum(sHit(mask))/sum(mask))*100, (sum(nHit(mask))/sum(mask))*100];
        else
            plotData{3,c} = [0, 0];
        end
    end
    
    % --- 3. CALCULATE ROW MARGINALS (Right Column) ---
    for r = 1:2
        fVal = 2 - r;
        mask = (fixCat == fVal); % Ignores saccade category
        if any(mask)
            plotData{r,3} = [(sum(sHit(mask))/sum(mask))*100, (sum(nHit(mask))/sum(mask))*100];
        else
            plotData{r,3} = [0, 0];
        end
    end

    % --- 4. PLOT THE 3x3 GRID ---
    figure;
    
    rowLabels = {'Fix: alvo', 'Fix: distrator', ''};
    colLabels = {'S: alvo', 'S: distrator', ''};
    
    for r = 1:3
        for c = 1:3
            % Skip the bottom-right corner (3,3) as requested
            if r == 3 && c == 3, continue; end
            
            % Calculate linear subplot index in a 3x3 grid
            subplotIdx = (r - 1) * 3 + c;
            subplot(3, 3, subplotIdx);
            
            % Render the 2 bars
            b = bar(plotData{r,c}, 'FaceColor', 'flat');
            
            % Style individual bars: sHit vs nHit
            b.CData(1,:) = [0.2 0.6 0.8]; % Blue for sHit
            b.CData(2,:) = [0.8 0.4 0.2]; % Orange/Red for nHit
            
            % UI Settings
            set(gca, 'TickDir', 'out', 'Box', 'off', 'XTickLabel', {'', ''});
            ylim([0 100]);
            
            % Add labels conditionally for cleanliness
            if c == 1, ylabel({['\bf{', rowLabels{r}, '}\rm'], 'Acertos (%)'}, 'Interpreter', 'tex'); end
            if r == 1, title(colLabels{c}); end
%             if r == 3 || (r == 2 && c == 3), xlabel('Tipo de Hit'); end
            
            grid on;
        end
    end
    
    % Add a simple floating legend on the empty bottom-right subplot area
    subplot(3,3,9); axis off;
    hold on;
    patch([0 1 1 0], [0 0 1 1], [0.2 0.6 0.8], 'EdgeColor', 'none');
    text(1.5, 0.5, 'S', 'VerticalAlignment', 'middle', 'FontSize', 11);
    patch([0 1 1 0], [-1.5 -1.5 -0.5 -0.5], [0.8 0.4 0.2], 'EdgeColor', 'none');
    text(1.5, -1, 'N', 'VerticalAlignment', 'middle', 'FontSize', 11);
    ylim([-2 2]); xlim([0 5]);

    sgtitle('Efeito de categoria dos estímulos')
end
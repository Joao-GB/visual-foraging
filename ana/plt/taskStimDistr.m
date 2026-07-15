function taskStimDistr(trl, mat, distOrder)
    if isfield(mat, 'drP') && isfield(mat.drP, 'blue')
        baseColor = mat.drP.blue;
    elseif isfield(mat, 'drP') && isfield(mat.drP, 'darkBlue')
        baseColor = mat.drP.darkBlue;
    else
        baseColor = [0, 0.4470, 0.7410];
    end

    % Cada célula guarda as diferenças para o k-ésimo vizinho
    relX = cell(1, distOrder);
    relY = cell(1, distOrder);
    
    for i = 1:length(trl)
        pos = trl(i).stmPosPix;
        M = size(pos, 2);
        
        if M < distOrder + 1
            continue;
        end
        
        dX = pos(1,:)' - pos(1,:); 
        dY = pos(2,:)' - pos(2,:); 
        distMat = sqrt(dX.^2 + dY.^2);
        
        for j = 1:M
            % Ordena as distâncias
            [~, sortIdx] = sort(distMat(j, :));
            
            % O k-ésimo vizinho está em k+1
            for k = 1:distOrder
                tgt = sortIdx(k + 1);
                relX{k}(end+1) = dX(tgt, j);
                relY{k}(end+1) = dY(tgt, j);
            end
        end
    end
    
    cMap = [linspace(1, baseColor(1), 256)', ...
            linspace(1, baseColor(2), 256)', ...
            linspace(1, baseColor(3), 256)'];
            
    % Uma figura para cada ordem de distância
    for k = 1:distOrder
        x = relX{k};
        y = relY{k};
        
        if isempty(x), continue; end
        
        % Converte para coordenadas polares para os gráficos marginais
        r = sqrt(x.^2 + y.^2);
        theta = atan2(y, x);
        
        figName = sprintf('Distribuição Espacial - %dº Vizinho', k);
        figure('Name', figName, 'Color', 'w', 'Position', [100 + k*30, 100 + k*30, 800, 750]);
        
        t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
        
        % Mapa de Densidade 2D nas 2 colunas superiores
        ax1 = nexttile([1 2]);

        histogram2(ax1, x, y, 50, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'off', 'EdgeColor', 'none');
        
        colormap(ax1, cMap);
        cb = colorbar(ax1);
        ylabel(cb, 'Ocorrências', 'FontWeight', 'bold');
        
        axis(ax1, 'equal');
        set(ax1, 'TickDir', 'out', 'Box', 'off');
        xlabel(ax1, '\Delta X (pixels)');
        ylabel(ax1, '\Delta Y (pixels)');
        title(ax1, sprintf('Densidade Relativa: %dº Vizinho Mais Próximo', k), 'FontSize', 14, 'FontWeight', 'bold');
        grid(ax1, 'on');
        
        % Distribuição radial na esquerda embaixo
        ax2 = nexttile;
        histogram(ax2, r, 40, 'FaceColor', baseColor, 'EdgeColor', 'k', 'LineWidth', 1.2);
        set(ax2, 'TickDir', 'out', 'Box', 'off');
        xlabel(ax2, 'Distância Radial (pixels)');
        ylabel(ax2, 'Frequência');
        title(ax2, 'Marginal Radial', 'FontSize', 12, 'FontWeight', 'bold');
        
        % Distribuição angular
        ax3 = polaraxes(t); 
        ax3.Layout.Tile = 4;
        
        polarhistogram(ax3, theta, 36, 'FaceColor', baseColor, 'EdgeColor', 'k', 'LineWidth', 1.2);
        title(ax3, 'Marginal Angular', 'FontSize', 12, 'FontWeight', 'bold');
    end
end
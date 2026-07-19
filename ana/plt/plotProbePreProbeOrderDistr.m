function plotProbePreProbeOrderDistr(trl, mat)
    % 1. Extração de variáveis base
    probeIdx      = [trl.probeIdx];
    preProbeIdx   = [trl.preProbeIdx];
    nSaccProbeIdx = [trl.nSaccProbeIdx];
    
    nStims  = mat.tkP.nStims;
    nTrials = numel(probeIdx);
    
    % 2. Identificação de Cores (com fallback caso não existam)
    if isfield(mat, 'drP') && isfield(mat.drP, 'green')
        cGreen  = mat.drP.green;
        cOrange = mat.drP.orange;
    else
        cGreen  = [0.4660 0.6740 0.1880];
        cOrange = [0.8500 0.3250 0.0980];
    end
    
    % 3. Construção da matriz lógica de histórico [nTrials x nStims]
    forHistIdx = {trl.forHistIdx};
%     histLen    = {trl.forHistLen};
    histBool   = repmat({zeros(1, nStims)}, nTrials, 1);
    
    % Usamos UniformOutput=false e depois convertemos para matriz
%     histIdxList = cellfun(@(x,y,z) assignOnes(x,y,z), forHistIdx(:), histLen(:), histBool(:), 'UniformOutput', false);
    histIdxList = cellfun(@(x,z) assignOnes(x,z), forHistIdx(:), histBool(:), 'UniformOutput', false);
    histIdx = cell2mat(histIdxList); % Matriz lógica N x nStims (1 = no histórico)

    % Vetores para armazenar os ranks resultantes
    rankSaccAll  = nan(1, nTrials);
    rankSaccVal  = nan(1, nTrials);
    rankNSaccAll = nan(1, nTrials);
    rankNSaccVal = nan(1, nTrials);
    
    % 4. Loop principal para calcular o ranking de distâncias
    for i = 1:nTrials
        pos = trl(i).stmPosPix; % 2 x nStims
        pp  = preProbeIdx(i);
        p   = probeIdx(i);
        ns  = nSaccProbeIdx(i);
        
        % Distância Euclidiana do pre-probe para TODOS os outros estímulos
        dX = pos(1,:) - pos(1,pp);
        dY = pos(2,:) - pos(2,pp);
        dists = sqrt(dX.^2 + dY.^2);
        
        %% --- FIGURA 1: SACCADIC PROBE ---
        % Subplot 1: Considera todos os estímulos (Rank = quantos são estritamente mais próximos + 1)
        rankSaccAll(i) = sum(dists < dists(p)) + 1;
        
        % Subplot 2: Exclui o histórico (apenas posições ~histIdx)
        validSacc = ~histIdx(i,:); 
        rankSaccVal(i) = sum(dists(validSacc) < dists(p)) + 1;
        
        %% --- FIGURA 2: NON-SACCADIC PROBE ---
        % Subplot 1: Considera todos os estímulos
        rankNSaccAll(i) = sum(dists < dists(ns)) + 1;
        
        % Subplot 2: Exclui o histórico E exclui o probe sacádico
        validNSacc = ~histIdx(i,:);
        validNSacc(p) = 0; 
        rankNSaccVal(i) = sum(dists(validNSacc) < dists(ns)) + 1;
    end
    
    % 5. Renderização
    renderRankFig('Ranking de distâncias entre pré-probe e probe sacádico', 'Pre-probe vs sacc probe dist ranking', ...
        rankSaccAll, rankSaccVal, ...
        'Todos estímulos', 'Sem estímulos do histórico', cGreen);
        
    renderRankFig('Ranking de distâncias entre pré-probe e probe não-sacádico', 'Pre-probe vs non-sacc probe dist ranking', ...
        rankNSaccAll, rankNSaccVal, ...
        'Todos estímulos', 'Sem estímulos do histórico e sacádico', cOrange);
end

%% --- HELPER FUNCTIONS ---

function z = assignOnes(x, z)
    % Verifica se y > 0 para evitar erros de indexação (x(1:0))
    if ~isempty(x) %&& y > 0
%         z(x(1:y)) = 1;
        z(x) = 1;
    end
end

function renderRankFig(figTitle, plotTitle, data1, data2, sub1Title, sub2Title, barColor)
    figure('Name', figTitle, 'Color', 'w', 'Position', [150 150 900 400]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    % Subplot 1
    ax1 = nexttile;
    % BinMethod 'integers' garante que cada barra represente exatamente 1 número inteiro de rank
    histogram(ax1, data1, 'BinMethod', 'integers', 'FaceColor', barColor, 'EdgeColor', 'w', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
    set(ax1, 'TickDir', 'out', 'Box', 'off');
    xlabel('Rank de Distância'); ylabel('Frequência (Trials)');
    title(sub1Title, 'FontSize', 12);
    
    % Subplot 2
    ax2 = nexttile;
    histogram(ax2, data2, 'BinMethod', 'integers', 'FaceColor', barColor, 'EdgeColor', 'w', 'LineWidth', 1.2, 'FaceAlpha', 0.9);
    set(ax2, 'TickDir', 'out', 'Box', 'off');
    xlabel('Rank de Distância');
    title(sub2Title, 'FontSize', 12);
    
    % Unifica limite X para facilitar a comparação visual
    maxRank = max([data1, data2]);
    if ~isnan(maxRank) && maxRank > 0
        xlim(ax1, [0.5, maxRank + 0.5]);
        xlim(ax2, [0.5, maxRank + 0.5]);
    end
    
    title(t, plotTitle, 'FontSize', 14, 'FontWeight', 'bold');
end




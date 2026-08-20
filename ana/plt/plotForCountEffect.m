function plotForCountEffect(trl, ~)
% Considera apenas os trials em que não há repetições/revisitas
    hasNoRep    = ~[trl.forHistHasRep];
    histLen     = {trl.forHistLen};      histLen = histLen(hasNoRep);
    ROILimsTime = {trl.ROILimsTime}; ROILimsTime = ROILimsTime(hasNoRep);
    stmLimsTime = {trl.stmLimsTime}; stmLimsTime = stmLimsTime(hasNoRep);
    
    % m.q. forHistFixDur, só que restrito aos índices de HistLen
    counts = cellfun(@(x,y,z)[diff([z(:,1) x(:,2:y)]); 1:y],ROILimsTime,histLen, stmLimsTime,'UniformOutput',false);
    counts = [counts{:}];

    histLen = [histLen{:}];
    
    forHit = {trl.forProbeHit}; forHit = forHit(hasNoRep);
    forHit = [forHit{:}];
    
    % Cria a figura e o layout
    figure('Name', 'Foraging count effects', 'Color', 'w', 'Position', [100, 100, 1000, 450]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    % Boxplots de duração da fixação por índice
    ax1 = nexttile;
    
    % Filtra NaNs para evitar erros no boxplot
    validCounts = ~isnan(counts(1,:)) & ~isnan(counts(2,:));
    cleanCounts = counts(:, validCounts);
    
    % O boxplot agrupa automaticamente a linha 1 (duração) pelos valores da linha 2 (índice)
    boxplot(cleanCounts(1,:), cleanCounts(2,:), ...
        'Colors', 'k', ...       % Caixas pretas para manter consistência
        'Symbol', 'k.');         % Outliers como pontos pretos
    
    % Estética do eixo
    set(ax1, 'TickDir', 'out', 'Box', 'off');
    xlabel('Índice da fixação na tentativa');
    ylabel('Duração (s)');
    title('Duração por índice da fixação');
    
    
    % Acurácia por histLen
    ax2 = nexttile;
    
    % Filtra NaNs
    validHist = ~isnan(histLen) & ~isnan(forHit);
    hData = histLen(validHist);
    fData = forHit(validHist);
    
    % Encontra os tamanhos de histórico únicos e calcula a acurácia para cada
    uHist = unique(hData);
    acc = zeros(size(uHist));
    nVals = zeros(size(uHist));
    
    for i = 1:length(uHist)
        idx = (hData == uHist(i));
        acc(i) = mean(fData(idx)) * 100; % Convertendo acurácia para % (0-100)
        nVals(i) = sum(idx);             % Quantidade de dados
    end
    
    % Plota as barras
    b = bar(uHist, acc, 0.6);
    
    % Aplica o estilo de forrageamento (barras brancas, bordas pretas grossas)
    b.FaceColor = 'w';
    b.EdgeColor = 'k';
    b.LineWidth = 1.5;
    
    % Estética do eixo
    set(ax2, 'TickDir', 'out', 'Box', 'off');
    ylim([0 115]); % Deixa 15% de respiro no topo para não cortar os textos "n="
    
    % Garante que o eixo X mostre apenas inteiros (caso pule algum valor de histLen)
    xticks(uHist);
    
    xlabel('Estímulos forrageados');
    ylabel('Acurácia (%)');
    title('Acertos por tamanho do histórico de forrageamento');
    
    % Adiciona os textos indicando o "n" acima de cada barra
    for i = 1:length(uHist)
        text(uHist(i), acc(i), sprintf('n=%d', nVals(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'Color', [0.35 0.35 0.35]);
    end
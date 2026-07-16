function plotForRecencyEffect(trl, ~)
    forHit     = [trl.forProbeHit];
    forRecency = [trl.forProbeRecency];

    uRec = unique(forRecency);

    acc   = zeros(size(uRec));
    nVals = zeros(size(uRec));

    for i = 1:length(uRec)
        idx = (forRecency == uRec(i));
        acc(i) = mean(forHit(idx)) * 100;
        nVals(i) = sum(idx);
    end

    figure('Name', 'Foraging recency effect', 'Color', 'w', 'Position', [150, 150, 700, 450]);

    b = bar(uRec, acc, 0.6);

    b.FaceColor = 'w';
    b.EdgeColor = 'k';
    b.LineWidth = 1.5;

    set(gca, 'TickDir', 'out', 'Box', 'off');
    xticks(uRec);
    grid on;

    xlabel('Recência em relação ao pré-probe');
    ylabel('Acurácia (%)');
    title('Efeito de recência no Desempenho', 'FontSize', 13, 'FontWeight', 'bold');

    ylim([0 115]); 

    for i = 1:length(uRec)
        if nVals(i) > 0
            text(uRec(i), acc(i), sprintf('n=%d', nVals(i)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'Color', [0.35 0.35 0.35]);
        end
    end

end
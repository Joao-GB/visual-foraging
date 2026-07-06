function plotPSAdurations(trl, drP)
    stimStart = arrayfun(@(s) s.phaseLimsTime(3,1), trl);
    stimEnd = arrayfun(@(s) s.phaseLimsTime(3,2), trl);
    
    fixStart = stimStart - arrayfun(@(s) s.P3FixDurPerPhase(2), trl);
    saccStart = stimEnd + arrayfun(@(s) s.P3FixDurPerPhase(4), trl);
    
    stimStart = stimStart(:); stimEnd = stimEnd(:);
    fixStart  = fixStart(:);   saccStart = saccStart(:);

    % Coluna 1: alinha ao início do estímulo (P3)
    c1_fix  = fixStart - stimStart;
    c1_end  = stimEnd - stimStart;
    
    % Coluna 2 e figura da direita: alinha ao início da sacada
    c2_fix  = fixStart - saccStart;
    c2_end  = stimEnd - saccStart;

    c3_stim = stimStart - saccStart;

    % Figura 2x3
    figure('Name', 'PSA Stim and Sacc Times', 'Color', 'w', 'Position', [50 50 1300 700]);
    t = tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'normal');
    
    colorLeft  = drP.orange;
    colorRight = drP.darkGreen;

    % Coluna 1: Encontra o maior limite absoluto para manter os zeros alinhados espacialmente
    maxC1 = ceil(max(max(c1_fix), max(c1_end))*10)/10; minC1 = floor(min(min(c1_fix), min(c1_end))*10)/10;
    xLimsC1 = [minC1 maxC1];

    binEdgesC1 = minC1 : 0.01 : maxC1;

    % Top Left: fixStart - stimStart
    nexttile(1);
    histogram(c1_fix, binEdgesC1, 'FaceColor', colorLeft, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    xlim(xLimsC1); xline(0, '--k', 'LineWidth', 1.2);
    xlabel('Tempo desde início do estímulo (s)'); ylabel('Tentativas');
    title('Início da fixação');
    grid on; fixAxisOffset();

    % Bottom Left: stimEnd - stimStart
    nexttile(4);
    histogram(c1_end, binEdgesC1, 'FaceColor', colorLeft, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    xlim(xLimsC1); xline(0, '--k', 'LineWidth', 1.2);
    xlabel('Tempo desde início do estímulo (s)'); ylabel('Tentativas');
    title('Fim do estímulo');
    grid on; fixAxisOffset();

    % Coluna 2
    maxC2 = ceil(max(max(c2_fix), max(c2_end))*10)/10; minC2 = floor(min(min(c2_fix), min(c2_end))*10)/10;
    xLimsC2 = [minC2 maxC2];

    binEdgesC2 = minC2 : 0.02 : maxC2;

    % Top Middle: fixStart - saccStart
    nexttile(2);
    histogram(c2_fix, binEdgesC2, 'FaceColor', colorRight, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    xlim(xLimsC2); xline(0, '--k', 'LineWidth', 1.2);
    xlabel('Tempo desde início da sacada (s)');
    title('Início da fixação');
    grid on; fixAxisOffset();

    % Bottom Middle: stimEnd - saccStart
    nexttile(5);
    histogram(c2_end, binEdgesC2, 'FaceColor', colorRight, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    xlim(xLimsC2); xline(0, '--k', 'LineWidth', 1.2);
    xlabel('Tempo desde início da sacada (s)');
    title('Fim do estímulo');
    grid on; fixAxisOffset();

    % Coluna 3
    maxC3 = ceil(max(c3_stim)*10)/10; minC3 = floor(min(c3_stim)*10)/10;
    binEdgesC3 = minC3 : 0.01 : maxC3;

    % Right Side: stimStart - saccStart (Saccade aligned)
    nexttile(3, [2, 1]); 
    histogram(c3_stim, binEdgesC3, 'FaceColor', colorRight, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
    xlim([minC3 maxC3]); xline(0, '--k', 'LineWidth', 1.2); 
    
    
    xlabel('Tempo desde início da sacada (s)'); ylabel('Tentativas');
    title('Início do estímulo');
    grid on; fixAxisOffset();

    % Global figure title
    title(t, 'Tempo entre ruído rosa e fixação', 'FontSize', 14, 'FontWeight', 'bold');
end


function fixAxisOffset()
    % Helper function to set formatting parameters and apply the 
    % custom vertical padding to keep zero bins from clashing with the X-axis line.
    set(gca, 'TickDir', 'out', 'Box', 'off');
    drawnow; % Forces generation of limits
    yl = ylim;
    ylim([-0.003 * diff(yl), yl(2)]);
end



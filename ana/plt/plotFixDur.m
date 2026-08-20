function goodStm = plotFixDur(trl, mat, plotMode)
if nargin < 3, plotMode = false; end

goodFix1 = arrayfun(@(s) diff(s.fixLimsTime(:, ...
    s.phaseLimsTime(3,2) > s.fixLimsTime(2,:) & ...
    s.phaseLimsTime(1,2) <= s.fixLimsTime(1,:))), ...
    trl, 'UniformOutput', false);
goodFix = [goodFix1{:}];

goodStm1 = arrayfun(@(s) diff(s.ROILimsTime(:, ...
    s.phaseLimsTime(3,2) > s.ROILimsTime(2,:) & ...
    s.phaseLimsTime(1,2) <= s.ROILimsTime(1,:))), ...
    trl, 'UniformOutput', false);
goodStm = [goodStm1{:}];


fixDursPerCountsOfFixPerStmPh2 = cell(1,10);
for i = 1:numel(trl)
    trlGoodStmPh2 = find(2 == trl(i).stmPerPhase);
    N = numel(trlGoodStmPh2);

    trlFixDurs = diff(trl(i).fixLimsTime);

    trlFixPerGoodStmPh2 = cell(1,N);

    for j = 1:N
        aux = sum(trl(i).fixPerStm,2) > 0;
        trl(i).fixPerStm = trl(i).fixPerStm(aux,:);

        k = trlGoodStmPh2(j);

        trlFixPerGoodStmPh2{j} = find(trl(i).fixPerStm(k,:));

        n = numel(trlFixPerGoodStmPh2{j});

        fixDursPerCountsOfFixPerStmPh2{n} = ...
            [fixDursPerCountsOfFixPerStmPh2{n}; ...
             trlFixDurs(trlFixPerGoodStmPh2{j})];
    end
end

fixDursPerCountsOfFixPerStmPh3 = cell(1,10);
for i = 1:numel(trl)
    trlGoodStmPh3 = find(3 == trl(i).stmPerPhase);
    N = numel(trlGoodStmPh3);

    trlFixDurs = diff(trl(i).fixLimsTime);

    trlFixPerGoodStmPh3 = cell(1,N);

    for j = 1:N
        aux = sum(trl(i).fixPerStm,2) > 0;
        trl(i).fixPerStm = trl(i).fixPerStm(aux,:);

        k = trlGoodStmPh3(j);

        trlFixPerGoodStmPh3{j} = find(trl(i).fixPerStm(k,:));

        n = numel(trlFixPerGoodStmPh3{j});

        fixDursPerCountsOfFixPerStmPh3{n} = ...
            [fixDursPerCountsOfFixPerStmPh3{n}; ...
             trlFixDurs(trlFixPerGoodStmPh3{j})];
    end
end


%% Figure setup (1x2 layout for side-by-side comparison)

figure('Name', 'PSA Valid Durations', ...
    'Color', 'w', ...
    'Position', [100 100 1000 450]);

t = tiledlayout(1, 2, ...
    'TileSpacing', 'loose', ...
    'Padding', 'normal');


% Assign colors from the provided 'mat' struct

colorFix = mat.drP.darkGreen;
colorStm = mat.drP.green;
lineColor = repmat(mat.drP.blackGrey, [1 3]);


%% Left Tile: Fixation Durations

maxF = ceil(max(goodFix)*10)/10;
minF = max(0, floor(min(goodFix)*10)/10);

binEdgesF = minF : 0.02 : maxF;

nexttile(1);

histogram(goodFix, binEdgesF, ...
    'FaceColor', colorFix, ...
    'EdgeColor', 'w', ...
    'FaceAlpha', 0.7);

xlim([minF maxF]);

xline(mean(goodFix), '--', ...
    'Color', lineColor, ...
    'LineWidth', 1.5);

xline(median(goodFix), '-', ...
    'Color', lineColor, ...
    'LineWidth', 1.5);

xlabel('Duração (s)');
ylabel('Quantidade');
title('Fixações individuais');

grid on;
fixAxisOffset();


%% Right Tile: Stimulus Durations

maxS = ceil(max(goodStm)*10)/10;
minS = max(0, floor(min(goodStm)*10)/10);

binEdgesS = minS : 0.02 : maxS;

nexttile(2);

histogram(goodStm, binEdgesS, ...
    'FaceColor', colorStm, ...
    'EdgeColor', 'w', ...
    'FaceAlpha', 0.7, ...
    'HandleVisibility', 'off');

xlim([minS maxS]);

xline(mean(goodStm), '--', ...
    'Color', lineColor, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Média');

xline(median(goodStm), '-', ...
    'Color', lineColor, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Mediana');

legend('Location', 'northeast', 'Box', 'off');

xlabel('Duração (s)');
ylabel('Quantidade');
title('Janela de fixações em estímulos');

grid on;
fixAxisOffset();


%% Global Formatting

title(t, 'Distribuição de durações', ...
    'FontSize', 14, ...
    'FontWeight', 'bold');


%% ================================================================
%  Plot distributions according to number of fixations per stimulus
%  ================================================================

if plotMode

    maxD = 0;

    for C = {fixDursPerCountsOfFixPerStmPh2, fixDursPerCountsOfFixPerStmPh3}
        fixDurs = C{1};
    
        for i = 1:numel(fixDurs)
            if ~isempty(fixDurs{i})
                d = fixDurs{i}(:);
                d = d(isfinite(d));
    
                if ~isempty(d)
                    maxD = max(maxD, max(d));
                end
            end
        end
    end
    maxD = ceil(maxD*10)/10;

    auxPlot(fixDursPerCountsOfFixPerStmPh2, colorFix, 'Fixações por janela em estímulo (fase 2)', maxD);
    auxPlot(fixDursPerCountsOfFixPerStmPh3, colorFix, 'Fixações por janela em estímulo (fase 3)', maxD);
end

end

function auxPlot(fixDurs, colorFix, plotTitle, globalMaxD)

    % Minimum number of observations required to plot a distribution
    minN = 5;

    M = numel(fixDurs);

    validCell = false(1,M);

    for i = 1:M

        thisData = fixDurs{i};

        if ~isempty(thisData) && size(thisData,1) >= minN
            validCell(i) = true;
        end

    end

    validIdx = find(validCell);
    M = numel(validIdx);

    % Nothing to plot
    if isempty(validIdx)
        warning(['No fixation-count groups contain at least %d ' ...
                 'stimuli. No additional figure was created.'], minN);
        return
    end

    allDurations = [];

    for k = 1:numel(validIdx)
        i = validIdx(k);
        allDurations = [allDurations; ...
            fixDurs{i}(:)];
    end

    % Remove invalid values just in case
    allDurations = allDurations(isfinite(allDurations));

    minD = max(0, floor(min(allDurations)*10)/10);
    maxD = globalMaxD;

    % Same 20 ms bin width as the first figure
    binEdges = minD : 0.02 : maxD;

    % Make sure the last edge contains the maximum value
    if binEdges(end) < maxD
        binEdges(end+1) = maxD;
    end

    figure('Name', 'Fixation Durations by Fixation Count', ...
        'Color', 'w', ...
        'Position', [100 50 900 1000]);

    t2 = tiledlayout(numel(validIdx), 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');


    % -------------------------------------------------------------
    % Color palette
    %
    % Use the existing green colors as the basis and interpolate
    % between lighter and darker versions.
    % -------------------------------------------------------------

    baseColor = colorFix;

    % Create a palette ranging from lighter to darker versions
    % of the same color.
    colorPalette = zeros(M,3);

    for i = 1:M

        frac = (i-1) / max(1,M-1);

        % Blend with white for the lower fixation counts
        colorPalette(i,:) = ...
            (1-frac)*min(1, baseColor + 0.35) + ...
            frac*baseColor;

    end


    % -------------------------------------------------------------
    % Plot each fixation-count group
    % -------------------------------------------------------------

    for k = 1:numel(validIdx)

        i = validIdx(k);

        thisData = fixDurs{i};

        ax = nexttile;

        hold(ax, 'on');

        % Plot each column separately.
        % Column j corresponds to the j-th fixation made to the
        % stimulus.
        for j = 1:i

            thisFix = thisData(:,j);
            thisFix = thisFix(isfinite(thisFix));

            if isempty(thisFix)
                continue
            end

            histogram(ax, thisFix, binEdges, ...
                'FaceColor', colorPalette(j,:), ...
                'EdgeColor', 'w', ...
                'FaceAlpha', 0.35, ...
                'EdgeAlpha', 0.25, ...
                'DisplayName', sprintf('Fixação %d',j));

        end

        xlim(ax, [minD maxD]);

        xlabel(ax, 'Duração (s)');
        ylabel(ax, 'Quantidade');

        title(ax, sprintf('%d fixações por estímulo (N = %d)', ...
            i, size(thisData,1)));

        grid(ax, 'on');

        fixAxisOffset();

        % Only show a legend when there is more than one fixation
        if i > 1
            legend(ax, ...
                'Location', 'northeast', ...
                'Box', 'off');
        end

        hold(ax, 'off');

    end


    % -------------------------------------------------------------
    % Global title
    % -------------------------------------------------------------

    title(t2, ...
        plotTitle, ...
        'FontSize', 14, ...
        'FontWeight', 'bold');
end


function fixAxisOffset()

% Helper function to set formatting parameters and apply the
% custom vertical padding to keep zero bins from clashing with
% the X-axis line.

set(gca, ...
    'TickDir', 'out', ...
    'Box', 'off');

drawnow;

% Forces generation of limits
yl = ylim;

ylim([-0.003 * diff(yl), yl(2)]);

end
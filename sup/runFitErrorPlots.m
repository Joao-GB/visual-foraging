function results = runFitErrorPlots( ...
    marginalizeParams, ...
    NumTrials, ...
    paramsGen, ...
    stairLevel, ...
    modelParams, ...
    NumSim)

    vary = [numel(NumTrials) > 1, cellfun(@numel, paramsGen) > 1];
    nVary = sum(vary);

    if nVary == 0
        error('At least one quantity must vary.');
    elseif nVary > 2
        error('At most two quantities may vary.');
    end

    values = [{NumTrials}, paramsGen];
    names  = {'Trials', '\alpha', '\beta', '\gamma', '\lambda'};

    varyIdx = find(vary);
    len = cellfun(@numel, values(varyIdx));

    if numel(varyIdx) == 1
        xIndex = varyIdx;
        curveIndex = [];
        xValues = values{xIndex};
        curveValues = [];
    else
        if len(1) >= len(2)
            xIndex = varyIdx(1);
            curveIndex = varyIdx(2);
        else
            xIndex = varyIdx(2);
            curveIndex = varyIdx(1);
        end

        xValues = values{xIndex};
        curveValues = values{curveIndex};
    end

    nX = numel(xValues);
    if isempty(curveIndex)
        nCurves = 1;
    else
        nCurves = numel(curveValues);
    end

    bias      = zeros(nCurves, nX);
    rmse      = zeros(nCurves, nX);
    sdError   = zeros(nCurves, nX);
    lower95   = zeros(nCurves, nX);
    upper95   = zeros(nCurves, nX);
    expectedX = zeros(nCurves, nX);

    for c = 1:nCurves
        for k = 1:nX

            currentTrials = NumTrials;
            currentParams = zeros(1, 4);

            for ii = 1:4
                if numel(paramsGen{ii}) == 1
                    currentParams(ii) = paramsGen{ii};
                end
            end

            if xIndex == 1
                currentTrials = xValues(k);
            else
                currentParams(xIndex - 1) = xValues(k);
            end

            if ~isempty(curveIndex)
                if curveIndex == 1
                    currentTrials = curveValues(c);
                else
                    currentParams(curveIndex - 1) = curveValues(c);
                end
            end

            resultsMC = runFitMonteCarlo( ...
                marginalizeParams, ...
                currentTrials, ...
                currentParams, ...
                stairLevel, ...
                modelParams, ...
                NumSim, ...
                false);

            bias(c, k)      = resultsMC.bias;
            rmse(c, k)      = resultsMC.rmse;
            sdError(c, k)   = std(resultsMC.error);
            lower95(c, k)   = prctile(resultsMC.error, 2.5);
            upper95(c, k)   = prctile(resultsMC.error, 97.5);
            expectedX(c, k) = resultsMC.expectedX;
        end
    end

    titleText = sprintf('%d simulations\nstairLevel = %.2f\n', NumSim, stairLevel);

    if isempty(curveIndex)
        titleText = sprintf('%s%s varying\n', titleText, names{xIndex});
    else
        titleText = sprintf('%s%s varying (x-axis)\n%s varying (curves)\n', ...
            titleText, names{xIndex}, names{curveIndex});
    end

    for ii = 1:4
        if numel(paramsGen{ii}) == 1
            titleText = sprintf('%s%s = %.3g\n', ...
                titleText, names{ii+1}, paramsGen{ii});
        end
    end

    titleText = sprintf('%smodelParams = [%.3g %.3g %.3g]', ...
        titleText, modelParams(1), modelParams(2), modelParams(3));

    if numel(NumTrials) == 1
        titleText = sprintf('%s\nTrials = %d', titleText, NumTrials);
    end

    figure('Color', 'w', 'Position', [100 100 700 700]);

    ax1 = subplot(2, 1, 1);
    hold(ax1, 'on');
    grid(ax1, 'on');
    box(ax1, 'on');

    colors = lines(nCurves);

    for c = 1:nCurves
        fill(ax1, ...
            [xValues fliplr(xValues)], ...
            [upper95(c, :) fliplr(lower95(c, :))], ...
            colors(c, :), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.20, ...
            'HandleVisibility', 'off');

        plot(ax1, xValues, bias(c, :), ...
            'Color', colors(c, :), ...
            'LineWidth', 2);
    end

    yline(ax1, 0, 'k--');
    ylabel(ax1, 'Bias');
    title(ax1, titleText, 'Interpreter', 'tex');

    if ~isempty(curveIndex)
        legendStrings = compose('%s = %.3g', names{curveIndex}, curveValues);
        lgd = legend(ax1, legendStrings, 'Location', 'best');
        lgd.Interpreter = 'tex';
    end

    ax2 = subplot(2, 1, 2);
    hold(ax2, 'on');
    grid(ax2, 'on');
    box(ax2, 'on');

    for c = 1:nCurves
        plot(ax2, xValues, rmse(c, :), ...
            'Color', colors(c, :), ...
            'LineWidth', 2);
    end

    xlabel(ax2, names{xIndex}, 'Interpreter', 'tex');
    ylabel(ax2, 'RMSE');

    if ~isempty(curveIndex)
        lgd = legend(ax2, legendStrings, 'Location', 'best');
        lgd.Interpreter = 'tex';
    end

    linkaxes([ax1 ax2], 'x');

    results.x          = xValues;
    results.xIndex     = xIndex;
    results.xName      = names{xIndex};
    results.curves     = curveValues;
    results.curveIndex  = curveIndex;

    if isempty(curveIndex)
        results.curveName = '';
    else
        results.curveName = names{curveIndex};
    end

    results.bias      = bias;
    results.rmse      = rmse;
    results.sdError   = sdError;
    results.lower95   = lower95;
    results.upper95   = upper95;
    results.expectedX = expectedX;
end
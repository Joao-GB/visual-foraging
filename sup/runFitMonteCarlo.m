function results = runFitMonteCarlo( ...
    marginalizeParams, NumTrials, paramsGen, stairLevel, modelParams, NumSim, doPlot)

    if nargin < 7 || isempty(doPlot)
        doPlot = true;
    end

    expectedX = PAL_CumulativeNormal( ...
        paramsGen, ...
        min(stairLevel, 1 - paramsGen(4) - 1e-4), ...
        'inverse');

    estX = zeros(NumSim, 1);

    for k = 1:NumSim
        [~, ~, estX(k)] = runFitSimulation( ...
            marginalizeParams, ...
            NumTrials, ...
            paramsGen, ...
            stairLevel, ...
            modelParams);
    end

    err = estX - expectedX;

    levelLabel = sprintf('%.0f', 100 * stairLevel);

    if doPlot
        figure;
        histogram(err, 30);
        xlabel(sprintf('Estimated x_{%s} - True x_{%s}', levelLabel, levelLabel));
        ylabel('Count');
        title(sprintf('%d simulations (%d trials each)', NumSim, NumTrials));
        grid on;

        figure;
        plot(estX, '.');
        hold on;
        yline(expectedX, 'r', 'LineWidth', 2);
        xlabel('Simulation');
        ylabel(sprintf('Estimated x_{%s}', levelLabel));
        title(sprintf('Distribution of x_{%s} estimates', levelLabel));
        grid on;
    end

    fprintf('True x%s = %.3f\n', levelLabel, expectedX);
    fprintf('Mean estimate = %.3f\n', mean(estX));
    fprintf('Bias = %.3f\n', mean(err));
    fprintf('SD = %.3f\n', std(estX));
    fprintf('RMSE = %.3f\n', sqrt(mean(err.^2)));

    results.estX      = estX;
    results.error     = err;
    results.expectedX = expectedX;
    results.bias      = mean(err);
    results.sd        = std(estX);
    results.rmse      = sqrt(mean(err.^2));
end
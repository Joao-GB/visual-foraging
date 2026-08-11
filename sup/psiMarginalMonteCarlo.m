function results = psiMarginalMonteCarlo(marginalizeParams, NumTrials, paramsGen, stairLevel, NumSim, doPlot, options)

arguments
    marginalizeParams = [4]
    NumTrials = 150
    paramsGen = [-50 0.15 0.5 0.03]
    stairLevel = 0.75
    NumSim = 100
    doPlot = true

    options.stimRange = [-80 -15]
    options.stimGrain = 10;
    options.stimScale = 'linear'

    options.plotGrain  = 51
    options.alphaRange = [-100 -10]
    options.betaRange  = [.05 1]
    options.betaScale  = 'log'

    options.priorAlphaRange = [-80 -20]
    options.priorAlphaGrain = 60
    options.priorAlphaMean  = -50
    options.priorAlphaStd   = 20
    
    options.priorBetaRange  = [1/(2^5) 2^3]
    options.priorBetaGrain  = 40
    options.priorBetaMean   = .5
    options.priorBetaStdOrders = 1

    options.priorLambdaGrain = 11;
    options.priorLambdaRange = [.01 .1];
    options.priorLambdaMean = .03;
    options.priorLambdaConc = 20;
end

plotOption = 'none';

% True stimulus intensity corresponding to the requested stair level
expectedX = PAL_CumulativeNormal( ...
    paramsGen, ...
    min(stairLevel, 1 - paramsGen(4) - 1e-4), ...
    'inverse');

estX = zeros(NumSim,1);

for k = 1:NumSim
    [~,~,estX(k)] = psiMarginalSimulation1( ...
    marginalizeParams, plotOption, NumTrials, ...
    paramsGen, stairLevel, ...
    'stimRange', options.stimRange, ...
    'stimGrain', options.stimGrain, ...
    'stimScale', options.stimScale, ...
    'plotGrain', options.plotGrain, ...
    'alphaRange', options.alphaRange, ...
    'betaRange', options.betaRange, ...
    'betaScale', options.betaScale, ...
    'priorAlphaRange', options.priorAlphaRange, ...
    'priorAlphaGrain', options.priorAlphaGrain, ...
    'priorAlphaMean', options.priorAlphaMean, ...
    'priorAlphaStd', options.priorAlphaStd, ...
    'priorBetaRange', options.priorBetaRange, ...
    'priorBetaGrain', options.priorBetaGrain, ...
    'priorBetaMean', options.priorBetaMean, ...
    'priorBetaStdOrders', options.priorBetaStdOrders, ...
    'priorLambdaGrain', options.priorLambdaGrain, ...
    'priorLambdaRange', options.priorLambdaRange, ...
    'priorLambdaMean', options.priorLambdaMean, ...
    'priorLambdaConc', options.priorLambdaConc);
end

error = estX - expectedX;

% Label corresponding to the requested threshold
levelLabel = sprintf('%.0f', 100*stairLevel);

%% Histogram of estimation error
if doPlot
    figure;
    histogram(error,30);

    xlabel(sprintf('Estimated x_{%s} - True x_{%s}', ...
        levelLabel, levelLabel));
    ylabel('Count');
    title(sprintf('%d simulations (%d trials each)', ...
        NumSim, NumTrials));
    grid on;

    %% Estimates across simulations
    figure;
    plot(estX,'.');
    hold on;
    yline(expectedX,'r','LineWidth',2);

    xlabel('Simulation');
    ylabel(sprintf('Estimated x_{%s}',levelLabel));
    title(sprintf('Distribution of x_{%s} estimates',levelLabel));
    grid on;
end

%% Summary statistics
fprintf('True x%s = %.3f\n',levelLabel,expectedX);
fprintf('Mean estimate = %.3f\n',mean(estX));
fprintf('Bias = %.3f\n',mean(error));
fprintf('SD = %.3f\n',std(estX));
fprintf('RMSE = %.3f\n',sqrt(mean(error.^2)));

%% Output
results.estX      = estX;
results.error     = error;
results.expectedX = expectedX;
results.bias      = mean(error);
results.sd        = std(estX);
results.rmse      = sqrt(mean(error.^2));

end

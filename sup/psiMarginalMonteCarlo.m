function results = psiMarginalMonteCarlo(marginalizeParams,NumTrials,paramsGen,stairLevel,NumSim, doPlot)
    if nargin < 6 || isempty(doPlot)
        doPlot = true;
    end

    plotOption = 'none';

    expectedX = PAL_CumulativeNormal(paramsGen,min(stairLevel,1-paramsGen(4)-1e-4),'inverse');

    estX = zeros(NumSim,1);

    for k = 1:NumSim
        [~,~,estX(k)] = psiMarginalSimulation(marginalizeParams,plotOption,NumTrials,paramsGen,stairLevel);
    end

    error = estX - expectedX;

    % Label corresponding to the requested threshold
    levelLabel = sprintf('%.0f',100*stairLevel);

    %% Histogram of estimation error
    if doPlot
        figure;
        histogram(error,30);
    
        xlabel(sprintf('Estimated x_{%s} - True x_{%s}',levelLabel,levelLabel));
        ylabel('Count');
        title(sprintf('%d simulations (%d trials each)',NumSim,NumTrials));
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
function plotStaircase(RF)

for b = 1:numel(RF)

    trialNum = 1:length(RF(b).x);
    
    presentedSigma = -RF(b).x;
    
    figure;
    hold on;
    
    % Trajectory line
    plot(trialNum, presentedSigma, ...
        'k-', 'LineWidth', 1.5);
    
    % Correct trials
    plot(trialNum(RF(b).response == 1), ...
         presentedSigma(RF(b).response == 1), ...
         'ko', ...
         'MarkerFaceColor', 'k', ...
         'MarkerSize', 7);
    
    % Incorrect trials
    plot(trialNum(RF(b).response == 0), ...
         presentedSigma(RF(b).response == 0), ...
         'ko', ...
         'MarkerFaceColor', 'w', ...
         'MarkerSize', 7);
    
    % Final threshold estimate
    yline(-RF(b).mean, '--k', 'LineWidth', 2);
    
    xlabel('Trial');
    ylabel('Sigma');
    title('AMRF Threshold Estimation');
    
    set(gca,'FontSize',14);
    
    ylim([5 80]);
    xlim([1 length(trialNum)]);
end
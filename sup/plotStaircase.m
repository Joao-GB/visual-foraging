function plotStaircase(RF, xmin, xmax)
    for b = 1:numel(RF)
        trialNum = 1:length(RF(b).x);
        presentedSigma = -RF(b).x;
        
        figure; hold on;
        
        % Trajetória
        plot(trialNum, presentedSigma, 'k-', 'LineWidth', 1.5);
        
        % Trials corretos
        plot(trialNum(RF(b).response == 1), presentedSigma(RF(b).response == 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        
        % Trials incorretos
        plot(trialNum(RF(b).response == 0), presentedSigma(RF(b).response == 0), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
        
        % Estimativa final do limiar
        yline(-RF(b).mean, '--k', 'LineWidth', 2);
        
        xlabel('Trial');
        ylabel('Sigma');
        title('AMRF Threshold Estimation');
        
        set(gca,'FontSize',14);
        
        ylim([xmin xmax]);
        xlim([1 length(trialNum)]);
    end
end
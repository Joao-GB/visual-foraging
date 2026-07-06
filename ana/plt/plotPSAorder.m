function order = plotPSAorder(trl, drP)

    [hits, counts, dPrimeData] = getPSAorder(trl); 
    
    order.hits = hits; 
    order.counts = counts;
    order.d = dPrimeData;

    titleOrder  = 'Efeito de ordem das respostas';
    xlabelOrder = {["F", "S", "N"], ["S", "F", "N"], ["S", "N", "F"], ["S", "N"];
                   ["F", "N", "S"], ["N", "F", "S"], ["N", "S", "F"], ["N", "S"]};

    % Figura de acurácia
    figNameAcc = 'PSA order effect - Accuracy';
    plotTitleAcc = [titleOrder, ' (Acurácia)'];
    ylabelOrderAcc = {'Acertos (%)', 'Acertos (%)'};
    isSensitivityPlot = false;
    renderPSAOrder(order.hits, xlabelOrder, ylabelOrderAcc, figNameAcc, plotTitleAcc, drP, isSensitivityPlot);

    % Figura de d-prime
    figNameSens = 'PSA order effect - Sensitivity';
    plotTitleSens = [titleOrder, ' (Sensibilidade)'];
    ylabelOrderSens = {'Sensibilidade (d'')', 'Sensibilidade (d'')'};
    isSensitivityPlot = true;
    renderPSAOrder(order.d, xlabelOrder, ylabelOrderSens, figNameSens, plotTitleSens, drP, isSensitivityPlot);
    
end
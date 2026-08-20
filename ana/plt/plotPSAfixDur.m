function plotPSAfixDur(trl, mat)
    saccInt = -[trl.saccInterval];

    saccMedian = quantile(saccInt,.5);
    saccQuartiles = quantile(saccInt,[0.25, 0.5, 0.75]);

    figTitleMedian    = 'PSA sacc interval median split';
    figTitleQuartiles = 'PSA sacc interval quartile split';
    plotTitle = 'Efeito de intervalo de sacada após estímulo';
    
    plotPSAsplit(saccMedian, saccInt, trl, mat, figTitleMedian, plotTitle);
    plotPSAsplit(saccQuartiles, saccInt, trl, mat, figTitleQuartiles, plotTitle);
end
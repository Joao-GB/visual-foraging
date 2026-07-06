function plotPSAfixDur(trl, mat)
    saccInt = -[trl.saccInt1];

    saccMedian = quantile(saccInt,.5);
    saccQuartiles = quantile(saccInt,[0.25, 0.5, 0.75]);

    figTitleMedian    = 'PSA sacc interval median split';
    figTitleQuartiles = 'PSA sacc interval quartile split';
    plotTitle = 'Efeito da latência de sacada após estímulo';
    
    plotPSAsplit(saccMedian, saccInt, trl, mat, figTitleMedian, plotTitle);
    plotPSAsplit(saccQuartiles, saccInt, trl, mat, figTitleQuartiles, plotTitle);
end
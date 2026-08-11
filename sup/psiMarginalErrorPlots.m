function results = psiMarginalErrorPlots( ...
    marginalizeParams,...
    NumTrials,...
    paramsGen,...
    stairLevel,...
    NumSim, options)
    arguments 
        marginalizeParams = [4] 
        NumTrials = 150 
        paramsGen = [-50 0.15 0.5 0.03] 
        stairLevel = 0.75 
        NumSim = 100 
    
        options.stimRange = [-80 -15] 
        options.stimGrain = 10; 
        options.stimScale = 'linear' 
        
        options.plotGrain = 51 
        options.alphaRange = [-100 -10] 
        options.betaRange = [.05 1] 
        options.betaScale = 'log' 
        
        options.priorAlphaRange = [-80 -20] 
        options.priorAlphaGrain = 60 
        options.priorAlphaMean = -50 
        options.priorAlphaStd = 20 
    
        options.priorBetaRange = [1/(2^5) 2^3] 
        options.priorBetaGrain = 40 
        options.priorBetaMean = .5 
        options.priorBetaStdOrders = 1 
    
        options.priorLambdaGrain = 11; 
        options.priorLambdaRange = [.01 .1]; 
        options.priorLambdaMean = .03; 
        options.priorLambdaConc = 20; 
end

%% Determine varying quantities

vary = [numel(NumTrials)>1, cellfun(@numel,paramsGen)>1];

nVary = sum(vary);

if nVary == 0
    error('At least one quantity must vary.');
elseif nVary > 2
    error('At most two quantities may vary.');
end

values = [{NumTrials}, paramsGen];
names  = {'Trials','\alpha','\beta','\gamma','\lambda'};

varyIdx = find(vary);
len = cellfun(@numel,values(varyIdx));

if numel(varyIdx)==1

    xIndex = varyIdx;
    curveIndex = [];

    xValues = values{xIndex};
    curveValues = [];

else

    % longest vector becomes x-axis
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

%% Allocate

bias     = zeros(nCurves,nX);
rmse     = zeros(nCurves,nX);
sdError  = zeros(nCurves,nX);
lower95  = zeros(nCurves,nX);
upper95  = zeros(nCurves,nX);

%% Monte-Carlo loops

for c = 1:nCurves

    for k = 1:nX

        currentTrials = NumTrials;
        currentParams = zeros(1,4);

        % Fixed parameters
        for ii = 1:4
            if numel(paramsGen{ii})==1
                currentParams(ii)=paramsGen{ii};
            end
        end

        % x-axis quantity
        if xIndex==1
            currentTrials = xValues(k);
        else
            currentParams(xIndex-1)=xValues(k);
        end

        % curve quantity
        if ~isempty(curveIndex)

            if curveIndex==1
                currentTrials = curveValues(c);
            else
                currentParams(curveIndex-1)=curveValues(c);
            end

        end

        resultsMC = psiMarginalMonteCarlo(marginalizeParams, ... 
            currentTrials, ... 
            currentParams, ... 
            stairLevel, ... 
            NumSim, ... 
            false, ... 
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

        bias(c,k)     = resultsMC.bias;
        rmse(c,k)     = resultsMC.rmse;
        sdError(c,k)  = std(resultsMC.error);
        lower95(c,k)  = prctile(resultsMC.error,2.5);
        upper95(c,k)  = prctile(resultsMC.error,97.5);

    end

end

%% Figure title

titleText = sprintf('%d simulations\nstairLevel = %.2f\n', ...
    NumSim, stairLevel);

if isempty(curveIndex)

    titleText = sprintf('%s%s varying\n', ...
        titleText,names{xIndex});

else

    titleText = sprintf('%s%s varying (x-axis)\n%s varying (curves)\n', ...
        titleText,...
        names{xIndex},...
        names{curveIndex});

end

fixedNames = {'\alpha','\beta','\gamma','\lambda'};

for ii=1:4

    if numel(paramsGen{ii})==1

        titleText = sprintf('%s%s = %.3g\n',...
            titleText,...
            fixedNames{ii},...
            paramsGen{ii});

    end

end

if numel(NumTrials)==1

    titleText = sprintf('%sTrials = %d',titleText,NumTrials);

end

%% Plot

figure('Color','w','Position',[100 100 700 700]);

ax1=subplot(2,1,1);
hold on
grid on
box on

colors = lines(nCurves);

for c=1:nCurves

    fill([xValues fliplr(xValues)],...
         [upper95(c,:) fliplr(lower95(c,:))],...
         colors(c,:),...
         'EdgeColor','none',...
         'FaceAlpha',0.20);

    plot(xValues,bias(c,:),...
         'Color',colors(c,:),...
         'LineWidth',2);

end

yline(0,'k--')

ylabel('Bias')
title(titleText,'Interpreter','tex')

if ~isempty(curveIndex)

    legendStrings = compose('%s = %.3g', ...
        names{curveIndex},curveValues);

    lgd=legend(ax1,legendStrings,'Location','best');
    lgd.Interpreter='tex';

end

ax2=subplot(2,1,2);
hold on
grid on
box on

for c=1:nCurves

    plot(xValues,rmse(c,:),...
     'Color',colors(c,:),...
     'LineWidth',2);

end

xlabel(names{xIndex},'Interpreter','tex')
ylabel('RMSE')

if ~isempty(curveIndex)

    lgd=legend(ax2,legendStrings,'Location','best');
    lgd.Interpreter='tex';

end

linkaxes([ax1 ax2],'x')

%% Output

results.x = xValues;
results.xIndex = xIndex;
results.xName = names{xIndex};

results.curves = curveValues;
results.curveIndex = curveIndex;

if isempty(curveIndex)
    results.curveName = '';
else
    results.curveName = names{curveIndex};
end

results.bias = bias;
results.rmse = rmse;
results.sdError = sdError;
results.lower95 = lower95;
results.upper95 = upper95;
end
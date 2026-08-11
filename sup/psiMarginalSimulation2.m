function [estParams, expectedX, estX] = psiMarginalSimulation2(marginalizeParams, plotOption, NumTrials, paramsGen, stairLevel, options)
% Versão em que a PF muda no meio da simulação: o alpha incrementa em 10
% unidades

    % Set default arguments if not provided
    arguments
        marginalizeParams = [4]
        plotOption = 'final'
        NumTrials = 150
        paramsGen = [-50 0.15 0.5 0.03]
        stairLevel = 0.75

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
        [priorLogBetaMean, priorLogBetaStd] = logBetaPrior(options.priorBetaMean, options.priorBetaStdOrders);
        stimRange = linspace(options.stimRange(1), options.stimRange(2), options.stimGrain);

    paramsGen2 = paramsGen;
    sum10 = 10*(2*(rand > 0.5) - 1)
    paramsGen2(1) = paramsGen2(1) + sum10;

    halfTrl = NumTrials/2;

    expectedX = PAL_CumulativeNormal(paramsGen2, min(stairLevel, 1-paramsGen(4)-.0001), 'inverse');

    % Seed random number generator
    if exist('RandStream.m','file')
        s = RandStream.create('mt19937ar','seed','shuffle'); 
        RandStream.setGlobalStream(s);
    end

    % Determine method flag for plotting titles/logic
    if isempty(marginalizeParams)
        method = 1; % Psi+
    else
        method = 2; % Psi-marginal
    end
    marginalize = marginalizeParams;

    % Parameter extraction
    gamma = paramsGen(3);
%     PF = @(params, x, varargin) PAL_CumulativeNormal_Level([params(1), 10.^params(2), params(3), params(4)], x, stairLevel, varargin{:});
%     PF = @(params, x, varargin) PAL_CumulativeNormal_Level(params, x, stairLevel, varargin{:});
    PF = @(params, x, varargin) PAL_CumulativeNormal(params, x, varargin{:});

     
    trackSuspend = [];
    suspend = 0;

    gammas = gamma;
    lambdas = linspace(options.priorLambdaRange(1),options.priorLambdaRange(2),options.priorLambdaGrain);
    
    if any(lambdas == 0)
        lambdas(lambdas == 0) = eps;
    end

    %%%%%%% Visual Setup (only if plotting)
    doPlot = ~strcmp(plotOption, 'none');
    if doPlot
        alphas = linspace(options.alphaRange(1), options.alphaRange(2), options.plotGrain);

        if strcmp(options.betaScale, 'log')
            logBetas = linspace(log10(options.betaRange(1)), log10(options.betaRange(2)), options.plotGrain);
        else
%             betas = linspace(options.betaRange(1), options.betaRange(2), options.plotGrain);
        end

        models = 1;
        [a, b, l] = ndgrid(alphas,logBetas,lambdas);
        a = permute(a, [2 1 3]);
        b = permute(b, [2 1 3]);
        l = permute(l, [2 1 3]);
        vizLUT = squeeze(PAL_AMPM_CreateLUT(alphas,logBetas,gammas,lambdas,models,stimRange,PF,false));
        vizLUT = permute(vizLUT,[2 1 3 4 5]);
        
        % Corrected visual prior
        vizprior = PAL_pdfNormal(a, options.priorAlphaMean, options.priorAlphaStd).*PAL_pdfNormal(b, priorLogBetaMean, priorLogBetaStd);
        posterior = vizprior;
        totals = zeros(size(stimRange));
        
        fh = figure('units','pixels','position',[100 50 1100 700]);
        axes
        set(gca,'units','pixels','position',[1 1 1100 700])
        set(gca,'xlim',[0 1],'ylim',[0 1])
        text(.1,.87,'(a)','fontsize',14)
        text(.48,.87,'(b)','fontsize',14)
        text(.66,.87,'(c)','fontsize',14)
        text(.84,.87,'(d)','fontsize',14)
        text(.46,.59,'(e)','fontsize',14)
        text(.04,.25,'(f)','fontsize',14)
        
        if method == 1
            text(.03,.97,'Psi+ method','fontsize',16)
            text(.03,.93,'Method maintains posterior in (a) and selects placement such as to minimize expected entropy in (a).','fontsize',12)
        elseif method == 2
            text(.03,.97,'Psi-marginal method','fontsize',16)
            text(.03,.93,'Method maintains posterior in (a) but selects placement such as to minimize expected entropy in (b).','fontsize',12)
        end
        axis off
        set(gca,'handlevisibility','off')
    end

    %%%%%%% Algorithm Setup
    priorAlphaRange = linspace(options.priorAlphaRange(1), options.priorAlphaRange(2), options.priorAlphaGrain);
    if strcmp(options.betaScale, 'log')
        priorLogBetaRange  = linspace(log10(options.priorBetaRange(1)), log10(options.priorBetaRange(2)), options.priorBetaGrain);
    else
%         priorBetaRange  = linspace(options.priorBetaRange(1), options.priorBetaRange(2), options.priorBetaGrain);
    end
    priorLambdaRange = lambdas;

    S = warning('QUERY', 'PALAMEDES:AMPM_setupPM:priorTranspose');
    warning('off','PALAMEDES:AMPM_setupPM:priorTranspose');

    PM = PAL_AMPM_setupPM('priorAlphaRange',single(priorAlphaRange),...
        'priorBetaRange',single(priorLogBetaRange),'priorGammaRange',single(gamma),...
        'priorLambdaRange',single(priorLambdaRange), 'numtrials',NumTrials, 'PF' , PF,...
        'stimRange',single(stimRange),'marginalize',marginalize);

    prior = PAL_pdfNormal(PM.priorAlphas, options.priorAlphaMean, options.priorAlphaStd);
    prior = prior.*PAL_pdfNormal(PM.priorBetas, priorLogBetaMean, priorLogBetaStd);
    prior = prior.*PAL_pdfBeta(PM.priorLambdas, options.priorLambdaMean, options.priorLambdaConc,'meanandconcentration'); 
    prior = prior./sum(sum(sum(sum(prior))));

    PM = PAL_AMPM_setupPM(PM,'prior',prior);
    warning(S);

    %%%%%%% Trial Loop
    count = 0;
    while PM.stop ~= 1
        
        % Simulate observer
        if count >= halfTrl
            response = rand(1) < PAL_CumulativeNormal(paramsGen2, PM.xCurrent);
        else
            response = rand(1) < PAL_CumulativeNormal(paramsGen, PM.xCurrent);
        end
        
        % Update visualization arrays if plotting is active
        if doPlot
            trackSuspend(length(PM.response)+1) = suspend; %#ok<AGROW> 
            idx = find(single(stimRange) == single(PM.xCurrent));
            totals(idx) = totals(idx) + 1;
            
            idxLUT = find(single(PM.stimRange) == single(PM.xCurrent));
            if response == 1
                posterior = posterior.*vizLUT(:,:,:,idxLUT);
            else
                posterior = posterior.*(1-vizLUT(:,:,:,idxLUT));
            end
            posterior = posterior./sum(sum(sum(posterior)));
        end
        
        % Update PM structure based on response
        PM = PAL_AMPM_updatePM(PM,response,'fixLapse',suspend);
        
        % Plot if 'all' is selected
        if strcmp(plotOption, 'all')
            renderPlots();
        end
        count = count + 1;
    end

    % Plot if 'final' is selected
    if strcmp(plotOption, 'final')
        renderPlots();
    end

    % Extract final estimates
    estParams = [ ...
    PM.threshold(end), ...
    10.^PM.slope(end), ...
    PM.guess(end), ...
    PM.lapse(end)];

    
    %%%%%%% Nested Plotting Function
    function renderPlots()
        drawnow
        clf(fh) % Clear specific figure
        
        % Re-draw labels
        axes('units','pixels','position',[1 1 1100 700])
        set(gca,'xlim',[0 1],'ylim',[0 1])
        text(.1,.87,'(a)','fontsize',14)
        text(.48,.87,'(b)','fontsize',14)
        text(.66,.87,'(c)','fontsize',14)
        text(.84,.87,'(d)','fontsize',14)
        text(.46,.59,'(e)','fontsize',14)
        text(.04,.25,'(f)','fontsize',14)
        if method == 1
            text(.03,.97,'Psi+ method','fontsize',16)
        elseif method == 2
            text(.03,.97,'Psi-marginal method','fontsize',16)
        end
        axis off
        set(gca,'handlevisibility','off')
        
        % (a): full 3-D posterior
        axes('units','pixels','position',[75 225 400 400]);
        slice(a,b,l,PAL_Scale0to1(posterior),[],[],lambdas);
        shading flat;
        if ~exist('OCTAVE_VERSION', 'builtin')
            alpha('color');
            am = linspace(.25,1,64);
            alphamap(am);
        end
        axis([min(alphas) max(alphas) min(logBetas) max(logBetas) min(lambdas) max(lambdas)])
        set(gca,'cameraposition',[-2.52 -1.57 0.1])
        set(gca,'xgrid', 'off','ygrid','off','zgrid','off')
        set(gca,'fontsize',8)
        xlabel('alpha','fontsize',12)
        ylabel('log_{10}(beta)','fontsize',12)
        zlabel('lambda','fontsize',12)
        
        % (b):
        if method == 1 || (method == 2 && sum(marginalize == 2))
            axes('units','pixels','position',[560 470 125 125]);
            toplot = squeeze(sum(sum(posterior(:,:,:),1),3))/max(squeeze(sum(sum(posterior(:,:,:),1),3)));
            plot(alphas,toplot,'k-')
            xlim_vals = [min(alphas) max(alphas)];
            set(gca,'xlim',xlim_vals)
            ylim_vals = [0 1.2];
            set(gca,'ylim',ylim_vals);
            xlabel('alpha','fontsize',10);
            text(xlim_vals(1),ylim_vals(2),'log_{10}(beta), \lambda marginalized','fontsize',10,'verticalalignment','bottom')
            set(gca,'xdir','normal')
            set(gca,'ydir','normal')  
        else
            axes('units','pixels','position',[560 470 125 125]);
            xlim_vals = [min(alphas) max(alphas)];
            set(gca,'xlim',xlim_vals)
            ylim_vals = [min(logBetas) max(logBetas)];
            set(gca,'ylim',ylim_vals);
            xlabel('alpha','fontsize',10);
            ylabel('log_{10}(beta)','fontsize',10)
            set(gca,'tickdir','out')        
            text(xlim_vals(1),ylim_vals(2),'\lambda marginalized','fontsize',10,'verticalalignment','bottom')
            toplot = (squeeze(sum(posterior,3))./max(max(sum(posterior,3))));
            axes('units','pixels','position',[560 470 125 125]);        
            image(alphas, logBetas, toplot,'cdatamapping','scaled')
            set(gca,'xdir','normal')
            set(gca,'ydir','normal')
            set(gca,'ytick',[], 'xtick',[],'fontsize',6)
        end
        
        % (c)
        if method == 1 || (method == 2 && sum(marginalize == 2))
            axes('units','pixels','position',[755 470 125 125]);
            toplot = squeeze(sum(sum(posterior(:,:,:),2),3))/max(squeeze(sum(sum(posterior(:,:,:),2),3)));
            plot(logBetas,toplot,'k-')
            xlim_vals = [min(logBetas) max(logBetas)];
            xlabel('log_{10}(beta)','fontsize',10);
            text(xlim_vals(2),1.2,'\alpha, \lambda marginalized','fontsize',10,'verticalalignment','bottom')
            set(gca,'xdir','reverse')
            set(gca,'xlim',xlim_vals)
            set(gca,'ylim',[0 1.2]);
        else
            axes('units','pixels','position',[755 470 125 125]);
            xlim_vals = [min(alphas) max(alphas)];
            set(gca,'xlim',xlim_vals)
            ylim_vals = [min(lambdas) max(lambdas)];
            set(gca,'ylim',ylim_vals);
            set(gca,'ytick',[0:.02:.1],'fontsize',8)
            xlabel('alpha','fontsize',10);
            ylabel('lambda','fontsize',10)
            set(gca,'tickdir','out')        
            text(xlim_vals(1),ylim_vals(2),'\beta marginalized','fontsize',10,'verticalalignment','bottom')
            toplot = (squeeze(sum(posterior,1))./max(max(sum(posterior,1))))';
            axes('units','pixels','position',[755 470 125 125]);        
            image(alphas, lambdas, toplot,'cdatamapping','scaled')
            set(gca,'xdir','normal')
            set(gca,'ydir','normal')
            set(gca,'ytick',[], 'xtick',[],'fontsize',6)
        end
        
        % (d)
        if method == 1 || (method == 2 && sum(marginalize == 2))
            axes('units','pixels','position',[950 470 125 125]);
            toplot = squeeze(sum(sum(posterior(:,:,:),1),2))/max(squeeze(sum(sum(posterior(:,:,:),1),2)));
            plot(lambdas,toplot,'k-')
            xlim_vals = [min(lambdas) max(lambdas)];
            set(gca,'xtick',xlim_vals,'ytick',[],'fontsize',8)
            xlabel('lambda','fontsize',10);
            text(xlim_vals(1),1.2,'\alpha, \beta marginalized','fontsize',10,'verticalalignment','bottom')
            set(gca,'xdir','normal')
            set(gca,'xlim',xlim_vals)
            set(gca,'ylim',[0 1.2]);
        else
            axes('units','pixels','position',[950 470 125 125]);
            xlim_vals = [min(logBetas) max(logBetas)];
            set(gca,'xlim',xlim_vals)
            ylim_vals = [0 .1];
            set(gca,'ylim',ylim_vals);
            set(gca,'ytick',[0:.02:.1],'fontsize',8)
            xlabel('log_{10}(beta)','fontsize',10);
            ylabel('lambda','fontsize',10)
            set(gca,'xdir','reverse')
            set(gca,'tickdir','out')        
            text(xlim_vals(2),ylim_vals(2),'\alpha marginalized','fontsize',10,'verticalalignment','bottom')
            toplot = (squeeze(sum(posterior,2))./max(max(sum(posterior,2))))';
            axes('units','pixels','position',[950 470 125 125]);        
            image(logBetas, lambdas, toplot,'cdatamapping','scaled')
            set(gca,'xdir','reverse')
            set(gca,'ydir','normal')
            set(gca,'ytick',[], 'xtick',[],'fontsize',6)
        end
        
        % (f)
        axes('units','pixels','position',[75 50 1000 125]);    
        t = 1:length(PM.x)-1;
        plot(t,PM.x(1:length(t)),'k');
        hold on
        plot(1:length(t),PM.threshold,'b-','LineWidth',2)
        plot(t(PM.response == 1 & trackSuspend == 0),PM.x(PM.response == 1 & trackSuspend == 0),'ko', 'MarkerFaceColor','k');
        plot(t(PM.response == 0 & trackSuspend == 0),PM.x(PM.response == 0 & trackSuspend == 0),'ko', 'MarkerFaceColor','w');
        
        for SR = 1:length(totals)
            if totals(SR) > 0
                plot(max(103,length(PM.x)+3),stimRange(SR),'ko','markerfacecolor','k','markersize',20*sqrt(totals(SR)./sum(totals)))
            end
        end
        
        minX = max(length(PM.x)-100, 0);
        maxX = max(length(PM.x)-100, 0)+105;
        minY = min(stimRange)-(max(stimRange)-min(stimRange))/10;
        maxY = max(stimRange)+(max(stimRange)-min(stimRange))/4;
        axis([minX maxX minY maxY])
        xlabel('Trial','fontsize',10);
        if strcmp(options.stimScale, 'log')
            ylabel('log(intensity) (\itx\rm)','fontsize',10);
        else
            ylabel('intensity (\itx\rm)','fontsize',10);
        end
        
        plot(minX + 5, maxY - (maxY-minY)/12, 'ko', 'markerfacecolor','k')
        text(double(minX + 6), double(maxY - (maxY-minY)/12), 'correct')    
        plot(minX + 12, maxY - (maxY-minY)/12, 'ko', 'markerfacecolor','w')
        text(double(minX + 13), double(maxY - (maxY-minY)/12), 'incorrect') 
           
        % (e)
        axes('units','pixels','position',[560 240 515 175]);
        hold on
        [SL, NP, OON] = PAL_PFML_GroupTrialsbyX(PM.x(1:length(PM.x)-1),PM.response,ones(size(PM.response)));
        for SR = 1:length(SL(OON~=0))
            plot(SL(SR),NP(SR)/OON(SR),'ko','markerfacecolor','k','markersize',20*sqrt(OON(SR)./sum(OON)))
        end
        axis([min(stimRange)-(max(stimRange)-min(stimRange))/10 max(stimRange)+(max(stimRange)-min(stimRange))/10 0 1]);
        
        % ML fit
        [~, MLindex] = PAL_findMax(PM.pdf);
        MLalpha = priorAlphaRange(MLindex(1));
        MLbeta = priorLogBetaRange(MLindex(2));
        MLlambda = priorLambdaRange(MLindex(4));
            
        if count >= halfTrl
            plot([min(stimRange):.01:max(stimRange)],PF(paramsGen,min(stimRange):.01:max(stimRange)),'g-','linewidth',2.5)
            plot([min(stimRange):.01:max(stimRange)],PF(paramsGen2,min(stimRange):.01:max(stimRange)),'k-','linewidth',2.5)
        else
            plot([min(stimRange):.01:max(stimRange)],PF(paramsGen,min(stimRange):.01:max(stimRange)),'k-','linewidth',2.5)
        end

        plot([min(stimRange):.01:max(stimRange)],PF([MLalpha 10.^MLbeta gamma MLlambda],min(stimRange):.01:max(stimRange)),'r-','linewidth',2)
        plot([min(stimRange):.01:max(stimRange)],PF([PM.threshold(end) 10.^PM.slope(end) gamma PM.lapse(end)],min(stimRange):.01:max(stimRange)),'b-','linewidth',2)
        
        if strcmp(options.stimScale, 'log')
            xlabel('log(intensity) (\itx\rm)','fontsize',12);
        else
            xlabel('intensity (\itx\rm)','fontsize',12);
        end
        ylabel('\psi(\itx\rm; \alpha, \beta, \gamma, \lambda)','fontsize',12);
        text(min(stimRange)+(max(stimRange)-min(stimRange))/4, .25,'Bayes','color','b','Fontsize',14)
        text(min(stimRange)+(max(stimRange)-min(stimRange))/4, .1,'ML (crude)','color','r','Fontsize',14)
        
        % Inset
        axes('units','pixels','position',[955 255 100 100]);    
        hold on
        plot([min(stimRange):.01:max(stimRange)],PF([0 1 0 0],min(stimRange):.01:max(stimRange)),'k-','linewidth',2) %#ok<*NBRAK2> 
        plot([min(stimRange):.01:max(stimRange)],PF([MLalpha 10.^MLbeta gamma MLlambda],min(stimRange):.01:max(stimRange)),'r-','linewidth',2)
        plot([min(stimRange):.01:max(stimRange)],PF([PM.threshold(end) 10.^PM.slope(end) gamma PM.lapse(end)],min(stimRange):.01:max(stimRange)),'b-','linewidth',2)
        set(gca,'xtick',[],'ytick',[0 1],'Fontsize',10);
        ylabel('F(\itx\rm)')
    end

    estX = PAL_CumulativeNormal(estParams, min(stairLevel, 1-estParams(4)-.0001), 'inverse');
end
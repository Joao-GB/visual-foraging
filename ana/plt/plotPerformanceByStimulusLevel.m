function plotPerformanceByStimulusLevel(X, cond, signal, response, varName, condNames, sgnNames, nbins, smoothWindow)
% Inputs (vetores de mesmo tamanho):
%   X         : preditor contínuo (e.g., força do estímulo)
%   cond      : índice da condição à qual o trial pertence
%   signal    : 1 = sinal presente; 0 = sinal ausente
%   response  : resposta do sujeito: 1 = afirma ter visto sinal
%
% Opcionais:
%   varName   : noma do preditor contínuo
%   condNames : nomes das condições
%   sgnNames  : nomes do sinal e do ruído
%   'nbins'        : number of X bins (default = 12)
%   'smoothWindow' : window for cumulative smoothing (default = 5)

arguments
    X
    cond
    signal
    response
    varName   = {'X', 'X (units)', 'extensive X name'}
    condNames = {'Cond1', 'Cond2'}
    sgnNames  = {'signal', 'noise'}
    nbins (1,1) double = 12
    smoothWindow (1,1) double = 5
end

cs = unique(cond);

% Verifica quais respostas estão corretas
correct = double(response == signal);

% Divide a variável contínua em bins
edges = linspace(min(X), max(X), nbins+1);
[~,~,bin] = histcounts(X, edges);

bin_centers = (edges(1:end-1) + edges(2:end))/2;


acc_inst  = nan(nbins,2,2); % [bin, cond, signal]
acc_cum   = nan(nbins,2,2);
hit_rate  = nan(nbins,2);
fa_rate   = nan(nbins,2);
dprime    = nan(nbins,2);
criterion = nan(nbins,2);

% Para cada condição e bin, ...
for c = 1:numel(cs)
    idxC = (cond == cs(c));

    for b = 1:nbins
        idx = (bin == b) & idxC;

        if sum(idx) < 3
            continue;
        end

        % ... encontra a acurácia instantânea, ...
        acc_inst(b,c,1) = mean(correct(idx & signal==1));
        acc_inst(b,c,2) = mean(correct(idx & signal==0));

        % ... a acumulada...
        idxCum = (bin <= b) & idxC;
        acc_cum(b,c,1) = mean(correct(idxCum & signal==1));
        acc_cum(b,c,2) = mean(correct(idxCum & signal==0));

        % ... e os valores de d' e c
        hits = sum(response(idx & signal==1)==1);
        misses = sum(response(idx & signal==1)==0);
        fas = sum(response(idx & signal==0)==1);
        crs = sum(response(idx & signal==0)==0);

        H = (hits + 0.5) / (hits + misses + 1);
        F = (fas + 0.5) / (fas + crs + 1);

        hit_rate(b,c) = H;
        fa_rate(b,c)  = F;

        dprime(b,c) = norminv(H) - norminv(F);
        criterion(b,c) = -0.5*(norminv(H) + norminv(F));
    end
end

% Suaviza um pouco os sinais acumulados
acc_cum = movmean(acc_cum, smoothWindow, 1, 'omitnan');

% Ajusta curvas sigmoides para cada subcondição
sigmoid = @(b, x) b(1) + (b(2) - b(1)) ./ (1 + exp(-(b(3) + b(4)*x)));
init = [0.5, 0.9, 0, 5];

fitParams = nan(2,2,4); % [cond, signal, 4 params]
sgn = [1 0];
for c = 1:numel(cs)
    for s = 1:2
        y = correct(cond == cs(c) & signal==sgn(s));
        x = X(cond == cs(c) & signal==sgn(s));
        model0 = fitnlm(x, y, sigmoid, init);
        fitParams(c,s,:) = model0.Coefficients.Estimate;
    end
end

% Plots
figure;

for c = 1:numel(cs)
    % Superior: acurácia
    subplot(2,numel(cs),c); hold on

    for s = 1:2
        idx = (cond == cs(c)) & (signal == sgn(s));
        xPts = X(idx);
        yPts = correct(idx); 
        if sgn(s) == 1
            plot(xPts, yPts, 'ko', 'MarkerSize', 6);
        else
            plot(xPts, yPts, 'ro', 'MarkerSize', 6);
        end
    end

    plot(bin_centers, acc_inst(:,c,1), '--k',  'LineWidth',1); % Sinal, instantânea
%     plot(bin_centers, acc_cum(:,c,1), '-k',  'LineWidth',3); % sinal, acumulada

    plot(bin_centers, acc_inst(:,c,2), '--r','LineWidth',1); % Ruído, instantânea
%     plot(bin_centers, acc_cum(:,c,2), '--r','LineWidth',3); % Ruído, acumulada

    xFine = linspace(min(X), max(X), 200);

    for s = 1:2
        b = squeeze(fitParams(c,s,:));
        if all(~isnan(b))
            yFit = sigmoid(b, xFine);
            if sgn(s) == 1
                plot(xFine, yFit, '-k', 'LineWidth',2);
            else
                plot(xFine, yFit, '-r', 'LineWidth',2);
            end
        end
    end

    xlabel(varName{2}); ylabel('Acurácia');
    title(sprintf('Condição %s: Acurácia em função da %s', condNames{c}, varName{3}));
    ylim([-0.05 1.05]);
    grid on;
    legend({sprintf('Acertos %s', sgnNames{1}),sprintf('Acertos %s', sgnNames{2}), ...
        sprintf('Inst. %s', sgnNames{1}),sprintf('Inst. %s', sgnNames{2}), ...
        sprintf('Ajuste %s',sgnNames{1}), sprintf('Ajuste %s',sgnNames{2})}, 'Location','best');

    % Inferior: d' e c
    subplot(2,numel(cs),c+numel(cs)); hold on

    plot(bin_centers, dprime(:,c), '-b','LineWidth',2);
    plot(bin_centers, criterion(:,c), '-m','LineWidth',2);
    
    xlabel(varName{2}); ylabel('Métricas de SDT');
    title(sprintf('Condição %s: ajustes psicométricos', condNames{c}));
    legend({'d''', 'c'}, 'Location', 'best');
    grid on;
end

end

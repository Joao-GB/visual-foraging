function [trlProps, results, eyeData, evTimes, goodTrlIdx] = foragingAnalysisOLD2(subj, ses)
    % Convenção:
    % vou chamar os estímulos de visto/fixado (F), pré-sacádico (pS) e 
    % não-sacádico (nS) e identificá-los como -1, 0 e 1 com respeito ao 
    % tempo em que são vistos quanto ao início da sacada

%     addpath(genpath('C:\Users\joaog\Documents\MATLAB\preCompiled_edfmex'));

    currFolder = fileparts(mfilename('fullpath')); parentFolder = fileparts(currFolder);
    addpath(genpath(fullfile(currFolder, 'dep'))); addpath(parentFolder);

    params = foragingParams;
    outFolder = fullfile(params.currFolder, params.outFolder);

    [mat, edf, sesStr] = foragingLoad(subj, ses, outFolder, params);

    mat.results.trialOrder(2,:,:) = mat.results.trialOrder(2,:,:) .* reshape(1:mat.tkP.nBlocks, 1, 1, []);

    %% Extrai as informações de cada trial, que já passaram por uma pré-seleção
    [trlProps, eyeData, evTimes] = foragingTrlProps(mat, edf, sesStr, subj);
    fprintf('\n----------------------------------\nNúmero de trials pré-selecionados: %d\n', sum(logical([trlProps.trlKeep])))

    %% Rejeita trials com base na quantidade de fixações durante o período de ruído rosa
    % e no período pós-ruído, devendo ser fixações distintas
    P3minDur = 0.02; P3maxDur = 0.180; maxOffset = .05;
    nTrl = numel(trlProps); 
    goodTrlIdx = zeros(1, nTrl);
    for i = 1:nTrl
        goodTrlIdx(i) = logical(trlProps(i).trlKeep) && trlProps(i).p3NFix >= 1 &&...
            numel(trlProps(i).p3p2FixIdx) == 1 && trlProps(i).pMOnlyNFix >= 1  &&...
            trlProps(i).p3p2FixIdx ~= trlProps(i).pMOnlyFixIdx(1) &&...
            trlProps(i).pinkNoiseDur >= P3minDur && trlProps(i).pinkNoiseDur <= P3maxDur && trlProps(i).p3FixOffsetDelay <= maxOffset;
    end

    psaTrlIdx = goodTrlIdx & [trlProps.p4StatusHas0and1];
    fprintf('Número final de trials selecionados: %d\n', sum(psaTrlIdx))

    psaTrl = trlProps(psaTrlIdx);


    %% Gráfico principal: Efeito pré-sacádico E tabela de contingência
    colNames = {'Trials com ruído seguido por sacada', 'Todos os trials'};
    rowNames = {'Proporção', 'Quantidade'};
    
    figure; x = [-1 0 1];

    fprintf('-----------------------\n TRIALS BONS\n')
    
    results.PSA.good = {}; results.PSA.all = {};
    [results.PSA.good.F.table, results.PSA.good.F.idx]   = contTable([psaTrl.p4AnsExpm1], [psaTrl.p4AnsObsm1]);
    results.PSA.good.F.correct = results.PSA.good.F.table(1,1)+results.PSA.good.F.table(2,2);
    results.PSA.good.F.total   = sum(results.PSA.good.F.table(:));

    fprintf('\n\nTabela: condição pré-sacádica (pS)')
    [results.PSA.good.pS.table, results.PSA.good.pS.idx, results.PSA.good.pS.d, results.PSA.good.pS.c] = ...
        contTable([psaTrl.p4AnsExp0], [psaTrl.p4AnsObs0], 1);
    results.PSA.good.pS.correct = trace(results.PSA.good.pS.table);
    results.PSA.good.pS.total   = sum(results.PSA.good.pS.table(:));


    fprintf('\n\nTabela: condição não sacádica (nS)')
    [results.PSA.good.nS.table, results.PSA.good.nS.idx, results.PSA.good.nS.d, results.PSA.good.nS.c] = ...
        contTable([psaTrl.p4AnsExp1], [psaTrl.p4AnsObs1], 1);
    results.PSA.good.nS.correct = trace(results.PSA.good.nS.table);
    results.PSA.good.nS.total   = sum(results.PSA.good.nS.table(:));
    

    counts = [results.PSA.good.F.correct results.PSA.good.pS.correct results.PSA.good.nS.correct; ...
        results.PSA.good.F.total results.PSA.good.pS.total results.PSA.good.nS.total];

    subplot(2,2,1);
    b = bar(x, counts(1,:)./counts(2,:), 'FaceColor', 'flat');
    b(1).CData(1,:) = [1 1 1];
    title(colNames{1}, 'Units','normalized', 'Position',[.5 1.00 0])
    yl1 = ylabel(rowNames{1}, 'FontWeight','bold');
    set(yl1, 'Units','normalized', 'Position',[-0.08 0.5 0]);
    set(gca, 'XTickLabel', []);

    subplot(2,2,3)
    b = bar(x, [counts(1,:);counts(2,:)-counts(1,:)], 'stacked', 'FaceColor', 'flat');
    b(1).CData(1,:) = [1 1 1];
    b(2).CData = repmat([.65 .65 .65], 3, 1);
    yl2 = ylabel(rowNames{2}, 'FontWeight','bold');
    set(yl2, 'Units','normalized', 'Position',[-0.08 0.5 0]);


    fprintf('-----------------------\n TODOS TRIALS\n')

    [results.PSA.all.F.table, results.PSA.all.F.idx]   = contTable([trlProps.p4AnsExpm1], [trlProps.p4AnsObsm1]);
    results.PSA.all.F.correct = results.PSA.all.F.table(1,1)+results.PSA.all.F.table(2,2);
    results.PSA.all.F.total   = sum(results.PSA.all.F.table(:));

    fprintf('\n\nTabela: condição pré-sacádica')
    [results.PSA.all.pS.table, results.PSA.all.pS.idx, results.PSA.all.pS.d, results.PSA.all.pS.c] = ...
        contTable([trlProps.p4AnsExp0], [trlProps.p4AnsObs0], 1);
    results.PSA.all.pS.correct = trace(results.PSA.all.pS.table);
    results.PSA.all.pS.total   = sum(results.PSA.all.pS.table(:));

    fprintf('\n\nTabela: condição não sacádica')
    [results.PSA.all.nS.table, results.PSA.all.nS.idx, results.PSA.all.nS.d, results.PSA.all.nS.c] = ...
        contTable([trlProps.p4AnsExp1], [trlProps.p4AnsObs1], 1);
    results.PSA.all.nS.correct = trace(results.PSA.all.nS.table);
    results.PSA.all.nS.total   = sum(results.PSA.all.nS.table(:));

    counts = [results.PSA.all.F.correct results.PSA.all.pS.correct results.PSA.all.nS.correct; ...
        results.PSA.all.F.total results.PSA.all.pS.total results.PSA.all.nS.total];

    subplot(2,2,2)
    b = bar(x, counts(1,:)./counts(2,:), 'FaceColor', 'flat');
    b(1).CData(1,:) = [1 1 1];
    title(colNames{2}, 'Units','normalized', 'Position',[.5 1.0 0])
    set(gca, 'XTickLabel', []);

    subplot(2,2,4)
    b = bar(x, [counts(1,:);counts(2,:)-counts(1,:)], 'stacked', 'FaceColor', 'flat');
    b(1).CData(1,:) = [1 1 1];
    b(2).CData = repmat([.65 .65 .65], 3, 1);

    sgtitle('Efeito pré-sacádico em tarefa de forrageamento');


    %% Efeito da ordem das perguntas
    % Não usei a infor do trlProps pois a função countFB já estava
    % estruturada para o trialFeedback, com 3 linhas por cell
    colNames = {'F primeiro', 'F no meio', 'F ausente ou ao fim', 'TOTAL'};
    rowNames = {'pS antes de nS', 'nS antes de pS'};
    figure;
    trlFB = reshape(mat.results.trialFeedback', 1, []);
    psaTrlFB = trlFB(psaTrlIdx);

    col1 = [1 1 1]; col2 = [.65 .65 .65];
        
    x = [-1 0 1];
    [fb_m101, N_m101, i_m101] = countFB(psaTrlFB, x, 3);
    foragingBarPlot(2,4,1, x, fb_m101, N_m101,[1 col1], col2);
    title(colNames{1}, 'Units','normalized', 'Position',[.5 1.00 0])
    ylabel(rowNames{1}, 'FontWeight','bold');
    

    x = [0 -1 1];
    [fb_0m11, N_0m11, i_0m11] = countFB(psaTrlFB, x, 3);
    foragingBarPlot(2,4,2, x, fb_0m11, N_0m11,[2 col1], col2);
    title(colNames{2}, 'Units','normalized', 'Position',[.5 1.00 0])

    x = [0 1];
    [fb_01, N_01, i_01] = countFB(psaTrlFB, x);
    foragingBarPlot(2,4,3, [x -1], fb_01, N_01,[3 col1], col2);
    title(colNames{3}, 'Units','normalized', 'Position',[.5 1.00 0])

%         i_01all = i_m101 | i_0m11 | i_01;
    fb_01all = [fb_m101(2)+fb_0m11(1)+fb_01(1) fb_m101(3)+fb_0m11(3)+fb_01(2)];
    N_01all  = [N_m101(2) + N_0m11(1)+ N_01(1)  N_m101(3)+ N_0m11(3)+ N_01(2)];
    foragingBarPlot(2,4,4, x, fb_01all, N_01all,[], col2);
    title(colNames{4}, 'Units','normalized', 'Position',[.5 1.00 0])


    x = [-1 1 0];
    [fb_m110, N_m110, i_m110] = countFB(psaTrlFB, x, 3);
    foragingBarPlot(2,4,5, x, fb_m110, N_m110,[1 col1], col2);
    ylabel(rowNames{2}, 'FontWeight','bold');

    x = [1 -1 0];
    [fb_1m10, N_1m10, i_1m10] = countFB(psaTrlFB, x, 3);
    foragingBarPlot(2,4,6, x, fb_1m10, N_1m10,[2 col1], col2);

    x = [1 0];
    [fb_10, N_10, i_10] = countFB(psaTrlFB, x); %#ok<*ASGLU> 
    foragingBarPlot(2,4,7, [x -1], fb_10, N_10,[3 col1], col2);

    fb_10all = [fb_m110(2)+fb_1m10(1)+fb_10(1) fb_m110(3)+fb_1m10(3)+fb_10(2)];
    N_10all  = [N_m110(2) + N_1m10(1)+ N_10(1)  N_m110(3)+ N_1m10(3)+ N_10(2)];
    foragingBarPlot(2,4,8, x, fb_10all, N_10all,[], col2);


    sgtitle('Efeito da ordem das perguntas na classificação dos estímulos');

    clear counts cPre cCurr cPost nPre nCurr nPost ...
        x b fb_10all fb_10 fb_m110 fb_1m10 fb_m101 fb_0m11 fb_01 fb_01all ...
        N_m101 N_0m11 N_01 N_01all N_m110 N_1m10 N_10 N_10all col1 col2 ...
        i_10 i_1m10 i_m110 i_01 i_0m11 i_m101

    %% Efeito da duração do ruído rosa: 
    % Esquerda: condição 1; Direita: condição 2
    % Cima: performance empírica, instantânea e ajustada; Baixo: medidas psicométricas
%     P3Dur = [psaTrl.pinkNoiseDur]; X = [P3Dur P3Dur];
%     sgn0 = [psaTrl.p4AnsExp0]; sgn1 = [psaTrl.p4AnsExp1]; sgn = [sgn0 sgn1];
%     rsp0 = [psaTrl.p4AnsObs0]; rsp1 = [psaTrl.p4AnsObs1]; rsp = [rsp0 rsp1];
%     cond0 = zeros(1,numel(sgn0)); cond1 = ones(1,numel(sgn1)); cond = [cond0 cond1];
%     varName = {'Duração ruído', 'Duração ruído (ms)', 'Duração do ruído rosa'};
%     condNames = {'pS', 'nS'};
%     sgnNames = {'alvo', 'distr.'};
% 
%     plotPerformanceByStimulusLevel(X, cond, sgn, rsp, varName, condNames, sgnNames);
%     sgtitle('Efeito da duração do ruído rosa no desempenho da tarefa de forrageamento')


    %% Efeito da distância entre o estímulo Curr e o Post perguntados
%     aux = [[psaTrl.p4PosDist1] [psaTrl.p3ToP4PosDist1]]; X = aux(1,:);
%     varName = {'Distância do estímulo', 'Distância (px)', 'Distância entre fixação e estímulo nS'};
%     condNames = {'pS a nS', 'último F a nS'};
%     plotPerformanceByStimulusLevel(X, [cond0 cond1], [sgn1 sgn1], [rsp1 rsp1], varName, condNames, sgnNames);
%     sgtitle('Efeito de distâncias no desempenho da tarefa de forrageamento')

    %% Heatmap das fixações ao redor do estímulo mais próximo quando há movimento ocular para fora do atual
%     fix.center = [psaTrl.pMFixPos]; fix.spread = [psaTrl.pMFixSpread];
% %     fixPointHeatmap([psaTrl.pMPos], fix, sgn0, rsp0, mat.txP.gabor.size_px, mat.txP.gabor.size_px*mat.prm.fixDistFactor3)
%     fixPointHeatmap([psaTrl.pMPos], fix.center, sgn0, rsp0, mat.txP.gabor.size_px, mat.txP.gabor.size_px*mat.prm.fixDistFactor3,'heatmap')
%         % É suficiente que saia do estímulo, e há dois casos: vai para
%         % outro estímulo ou não. Se não for, procuro qual o mais próximo,
%         % pois provavelmente havia intenção de visitá-lo
%         % Ou melhor, pego só os casos em que o olho foi para um estímulo,
%         % mas separo conforme acerots erros e alvo sem alvo


        %% Distribuição da duração das fixações, 3 curvas por plot (todas, alvo, não alvo)
        % em cima, direita considerando apenas as fixações antes do ruído
        % rosa, esquerda fixações durante PSA
        % (i.e., só P2)
        % embaixo, funde direita c esquerda e mantém os 3 plots

        %% Duração da fixação em relação ao número de estímulos já vistos


        % Distribuição da fixação em relação ao centro do estímulo, algo
        % que relacione a posição da fixação com a taxa de acerto


end

% Tabela de contingência:
% linhas:  sinal (presente, ausente)
% colunas: detecção (sim, não)
function [T,I, d, c] = contTable(E, O, showTable)
    if nargin < 3, showTable = 0; end 
    T = zeros(2,2); I = cell(2,2);
    % hit
    I{1,1} = E == 1 & O == 1;
    T(1,1) = sum(I{1,1});
    % false alarm
    I{1,2} = E == 0 & O == 1;
    T(1,2) = sum(I{1,2});
    % miss
    I{2,1} = E == 1 & O == 0;
    T(2,1) = sum(I{2,1});
    % rejeição correta
    I{2,2} = E == 0 & O == 0;
    T(2,2) = sum(I{2,2});

    zH = norminv(T(1,1)/(T(1,1)+T(2,1))); zFA = norminv(T(1,2)/(T(1,2)+T(2,2)));
    d = zH - zFA;
    c = -(zH + zFA)/2;

    if showTable
        fprintf('\n');
        fprintf('        Alvo presente     \n');
        fprintf('      -----------------     \n');
        fprintf('        Sim        Não       \n');
        fprintf('       -----      -----      \n');
        fprintf('Sim  %4d(H)   %4d(FA)  \n', T(1,1), T(1,2));
        fprintf('Não  %4d(M)   %4d(CR)  \n', T(2,1), T(2,2));
        fprintf('\n');
        
        % Totais
        fprintf('Totais: Alvo = %d, Distrator = %d\n', T(1,1)+T(2,1), T(1,2)+T(2,2));
        fprintf('d''  %.4f, c = %.4f\n', d, c);
    end
end
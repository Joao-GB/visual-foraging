function [trlProps, analysis, eyeData, evTimes] = foragingAnalysis(subj, ses)
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
    clear currFolder edf outFolder params parentFolder sesStr ses subj;

    trlProps = trlProps([trlProps.trlKeep]);

    fprintf('\n----------------------------------\nNúmero de trials pré-selecionados: %d\n', sum(logical([trlProps.trlKeep])))

    %% Rejeita trials com base no quanto de ruído rosa se vê
    P3SaccLatencyLims = -1*[mat.prm.pinkNoiseDur+.2 .020];
    P3SaccLatency = [trlProps.saccInt1]/1000;
    goodTrl = P3SaccLatency >= P3SaccLatencyLims(1) & P3SaccLatency <= P3SaccLatencyLims(2);

    trlProps = trlProps(goodTrl);

    fprintf('Número final de trials selecionados: %d\n', sum(goodTrl))



    %% Gráfico principal: Efeito pré-sacádico E tabela de contingência
    colNames = {'Trials válidos'};
    rowNames = {'Proporção', 'Quantidade'};
    
    figure; x = [-1 0 1];

    fprintf('-----------------------\n TRIALS BONS\n')
    
    analysis.PSA.good = {};
    [analysis.PSA.good.F.table, analysis.PSA.good.F.idx]   = contTable([trlProps.forProbeCat], [trlProps.forProbeResp]);
    analysis.PSA.good.F.correct = analysis.PSA.good.F.table(1,1)+analysis.PSA.good.F.table(2,2);
    analysis.PSA.good.F.total   = sum(analysis.PSA.good.F.table(:));

    fprintf('\n\nTabela: condição pré-sacádica (pS)')
    [analysis.PSA.good.pS.table, analysis.PSA.good.pS.idx, analysis.PSA.good.pS.d, analysis.PSA.good.pS.c] = ...
        contTable([trlProps.probeCat], [trlProps.probeResp], 1);
    analysis.PSA.good.pS.correct = trace(analysis.PSA.good.pS.table);
    analysis.PSA.good.pS.total   = sum(analysis.PSA.good.pS.table(:));


    fprintf('\n\nTabela: condição não sacádica (nS)')
    [analysis.PSA.good.nS.table, analysis.PSA.good.nS.idx, analysis.PSA.good.nS.d, analysis.PSA.good.nS.c] = ...
        contTable([trlProps.nSaccProbeCat], [trlProps.nSaccProbeResp], 1);
    analysis.PSA.good.nS.correct = trace(analysis.PSA.good.nS.table);
    analysis.PSA.good.nS.total   = sum(analysis.PSA.good.nS.table(:));
    

    counts = [analysis.PSA.good.F.correct analysis.PSA.good.pS.correct analysis.PSA.good.nS.correct; ...
        analysis.PSA.good.F.total analysis.PSA.good.pS.total analysis.PSA.good.nS.total];

    condName = ["Forrageado", "Sacádico", "Não-sacádico"];

        subplot(1,2,1);
        b = bar(counts(1,:)./counts(2,:), 'FaceColor', 'flat');
        b(1).CData(1,:) = [1 1 1];
        title(colNames{1}, 'Units','normalized', 'Position',[.5 1.00 0])
        yl1 = ylabel(rowNames{1}, 'FontWeight','bold');
        set(yl1, 'Units','normalized', 'Position',[-0.08 0.5 0]);
    %     set(gca, 'XTickLabel', []);
        set(gca, 'TickDir', 'out', 'Box', 'off')
        xticklabels(condName);

    subplot(1,2,2)
    b = bar(x, [counts(1,:);counts(2,:)-counts(1,:)], 'stacked', 'FaceColor', 'flat');
    b(1).CData(1,:) = [1 1 1];
    b(2).CData = repmat([.65 .65 .65], 3, 1);
    yl2 = ylabel(rowNames{2}, 'FontWeight','bold');
    set(yl2, 'Units','normalized', 'Position',[-0.08 0.5 0]);
    xticklabels(condName);

    sgtitle('Efeito pré-sacádico em tarefa de forrageamento');


    %% Efeito da ordem das perguntas
    % Não usei a infor do trlProps pois a função countFB já estava
    % estruturada para o trialFeedback, com 3 linhas por cell
    colNames = {'F primeiro', 'F no meio', 'F ausente ou ao fim', 'TOTAL'};
    rowNames = {'pS antes de nS', 'nS antes de pS'};
    figure;
    trlFB = reshape(mat.results.trialFeedback', 1, []);
    trlPropsFB = trlFB(goodTrl);

    col1 = [1 1 1]; col2 = [.65 .65 .65];
        
    x = [-1 0 1];
    [fb_m101, N_m101, i_m101] = countFB(trlPropsFB, x, 3);
    foragingBarPlot(2,4,1, x, fb_m101, N_m101,[1 col1], col2);
    title(colNames{1}, 'Units','normalized', 'Position',[.5 1.00 0])
    ylabel(rowNames{1}, 'FontWeight','bold');
    

    x = [0 -1 1];
    [fb_0m11, N_0m11, i_0m11] = countFB(trlPropsFB, x, 3);
    foragingBarPlot(2,4,2, x, fb_0m11, N_0m11,[2 col1], col2);
    title(colNames{2}, 'Units','normalized', 'Position',[.5 1.00 0])

    x = [0 1];
    [fb_01, N_01, i_01] = countFB(trlPropsFB, x);
    foragingBarPlot(2,4,3, [x -1], fb_01, N_01,[3 col1], col2);
    title(colNames{3}, 'Units','normalized', 'Position',[.5 1.00 0])

%         i_01all = i_m101 | i_0m11 | i_01;
    fb_01all = [fb_m101(2)+fb_0m11(1)+fb_01(1) fb_m101(3)+fb_0m11(3)+fb_01(2)];
    N_01all  = [N_m101(2) + N_0m11(1)+ N_01(1)  N_m101(3)+ N_0m11(3)+ N_01(2)];
    foragingBarPlot(2,4,4, x, fb_01all, N_01all,[], col2);
    title(colNames{4}, 'Units','normalized', 'Position',[.5 1.00 0])


    x = [-1 1 0];
    [fb_m110, N_m110, i_m110] = countFB(trlPropsFB, x, 3);
    foragingBarPlot(2,4,5, x, fb_m110, N_m110,[1 col1], col2);
    ylabel(rowNames{2}, 'FontWeight','bold');

    x = [1 -1 0];
    [fb_1m10, N_1m10, i_1m10] = countFB(trlPropsFB, x, 3);
    foragingBarPlot(2,4,6, x, fb_1m10, N_1m10,[2 col1], col2);

    x = [1 0];
    [fb_10, N_10, i_10] = countFB(trlPropsFB, x); %#ok<*ASGLU> 
    foragingBarPlot(2,4,7, [x -1], fb_10, N_10,[3 col1], col2);

    fb_10all = [fb_m110(2)+fb_1m10(1)+fb_10(1) fb_m110(3)+fb_1m10(3)+fb_10(2)];
    N_10all  = [N_m110(2) + N_1m10(1)+ N_10(1)  N_m110(3)+ N_1m10(3)+ N_10(2)];
    foragingBarPlot(2,4,8, x, fb_10all, N_10all,[], col2);


    sgtitle('Efeito da ordem das perguntas na classificação dos estímulos');



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
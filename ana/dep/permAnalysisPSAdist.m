function [p, T, permT] = permAnalysisPSAdist(allTrlProps, allTrlProbes, allSubj, nPerm, mat, mode, params, doPlot)
    if nargin < 5, doPlot = 0; end
    %% Numeração dos efeitos
    % ANOVA 2x2 com medidas repetidas em ambos os fatores:
    % 1: PSA;
    % 2: distância

    nSubj = numel(allSubj);

    %% --- Cria as máscaras de distância

    nCondFactor2 = 2;
    factor2 = cell(1, nSubj);
    for subj=1:nSubj
        [closeTrl, farTrl] = plotPSAStimTriangleProps1(allTrlProbes(subj).pre.pos, allTrlProbes(subj).sacc.pos, allTrlProbes(subj).nSacc.pos, mat.drP, mode,0);
        factor2{subj} = [closeTrl'; farTrl'];
    end

    dPerCond = zeros(nSubj,2,2);

    %% --- Obtém as estatísticas para os dados: 2 efeitos principais e uma interação
    % Repare que para cada sujeito pego os trials referentes a ele e filtro
    % conforme o fator 2. Não necessariamente a análise vai incluir todos
    % os trials relativos a cada sujeito. Mesmo procedimento para
    % permutação
    for i = 1:nSubj
        subj = allSubj(i);
        for j=1:nCondFactor2
            subjCondTrlProps = allTrlProps([allTrlProps.subjNum] == subj);
            subjCondTrlProps = subjCondTrlProps(factor2{i}(j,:));
            [main, ~] = getPSAeffect(subjCondTrlProps, 0);
            dPerCond(i,:, j) = [main.sacc.d, main.nSacc.d];
        end
    end
    Tfactor1 = sum(sum(dPerCond, [1 3]).^2);
    Tfactor2 = sum(sum(dPerCond, [1 2]).^2);
    % Para a interação eu subtraio ao longo do fator mais importante, mas
    % não deveria fazer diferença no resultado
    Tinter = sum(sum(diff(dPerCond, 1, 2), 1).^2);
    T = [Tfactor1 Tfactor2 Tinter];

    %% --- Realiza permutações e obtém a estatística para cada uma delas
    nStats = 3;
    permT = zeros(nStats,nPerm);
    
    for b = 1:nPerm
        %% (a) H_0,1: efeito 1 nulo. Embaralho condição sacádica, mantendo
        % separados os trials quanto aos demais fatores. Como só permuto dentro
        % de cada trial, não há porque separar de acordo com os outros fatores
        % antes de permutar
        permAllTrlProps = foragingPSApermute(allTrlProps);
        permDPerCond = zeros(nSubj, 2);
        for i = 1:nSubj
            subj = allSubj(i);
            for j=1:nCondFactor2
                permSubjCondTrlProps = permAllTrlProps([permAllTrlProps.subjNum] == subj);
                permSubjCondTrlProps =  permSubjCondTrlProps(factor2{i}(j,:));
                [main, ~] = getPSAeffect(permSubjCondTrlProps, 0);
                permDPerCond(i,:, j) = [main.sacc.d, main.nSacc.d];
            end
        end
        permT(1,b) = sum(sum(permDPerCond, [1 3]).^2);

        %% (b) H_0,2; efeito 2 nulo. Embaralho o close e o far

        permFactor2 =cell(1,nSubj);
        for i = 1:nSubj
            idx = find(any(factor2{i}, 1));
            labels = 1 + (factor2{i}(2, idx) == 1);
            labels = labels(randperm(numel(labels)));
            permFactor2{i} = factor2{i};
            permFactor2{i}(:, idx) = 0;
            permFactor2{i}(1, idx(labels == 1)) = 1;
            permFactor2{i}(2, idx(labels == 2)) = 1;
        end

        permDPerCond = zeros(nSubj, 2);
        for i = 1:nSubj
            subj = allSubj(i);
            for j=1:nCondFactor2
                permSubjCondTrlProps = allTrlProps([allTrlProps.subjNum] == subj);
                permSubjCondTrlProps = permSubjCondTrlProps(permFactor2{i}(j,:));
                [main, ~] = getPSAeffect(permSubjCondTrlProps, 0);
                permDPerCond(i,:, j) = [main.sacc.d, main.nSacc.d];
            end
        end
        permT(2,b) = sum(sum(permDPerCond, [1 2]).^2);

        %% (c) H_0 de interação: diz que a diferença dos níveis do fator 1 não depende do nível do fator 2
        % Parece que posso usar o mesmo esquema de permutação ao longo dos
        % labels do fator 2, o que muda é a estatística
        permT(3,b) = sum(sum(diff(permDPerCond, 1, 2), 1).^2);
    end

    p = zeros(1, nStats);
    for i = 1:nStats
        p(i) = (1 + sum(abs(permT(i,:)) >= abs(T(i)))) / (numel(permT(i,:)) + 1);
    end
end
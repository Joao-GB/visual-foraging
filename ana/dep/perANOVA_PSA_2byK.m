function [dPerCond, p, T, permT] = perANOVA_PSA_2byK(allTrlProps, subjFactor2, allSubj, nPerm, K)
    %% Numeração dos efeitos
    % ANOVA 2x2 com medidas repetidas em ambos os fatores:
    % 1: PSA;
    % 2: o outro efeito

    % subjFactor2 deve estar separado por sujeitos, com máscaras de índices
    % para cada nível

    nSubj = numel(allSubj);

    dPerCond = zeros(nSubj,2,K);

    %% --- Obtém as estatísticas para os dados: 2 efeitos principais e uma interação
    % Repare que para cada sujeito pego os trials referentes a ele e filtro
    % conforme o fator 2. Não necessariamente a análise vai incluir todos
    % os trials relativos a cada sujeito. Mesmo procedimento para
    % permutação
    for i = 1:nSubj
        subj = allSubj(i);
        for j=1:K
            subjCondTrlProps = allTrlProps([allTrlProps.subjNum] == subj);
            subjCondTrlProps = subjCondTrlProps(subjFactor2{i}(j,:));
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
            for j=1:K
                permSubjCondTrlProps = permAllTrlProps([permAllTrlProps.subjNum] == subj);
                permSubjCondTrlProps =  permSubjCondTrlProps(subjFactor2{i}(j,:));
                [main, ~] = getPSAeffect(permSubjCondTrlProps, 0);
                permDPerCond(i,:, j) = [main.sacc.d, main.nSacc.d];
            end
        end
        permT(1,b) = sum(sum(permDPerCond, [1 3]).^2);

        %% (b) H_0,2; efeito 2 nulo. Embaralho o close e o far

        permFactor2 =cell(1,nSubj);
        for i = 1:nSubj
            idx = find(any(subjFactor2{i}, 1));
            labels = ones(size(idx));
            for j=2:K
                labels = labels + (j-1)*(subjFactor2{i}(j, idx) == 1);
            end
            labels = labels(randperm(numel(labels)));
            permFactor2{i} = subjFactor2{i};
            permFactor2{i}(:, idx) = 0;
            for j=1:K
                permFactor2{i}(j, idx(labels == j)) = 1;
            end
        end

        permDPerCond = zeros(nSubj, 2);
        for i = 1:nSubj
            subj = allSubj(i);
            for j=1:K
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
function [pPSA] = foragingAnalysis(allSubj, searchFolder, nPerm)
    currFolder = fileparts(mfilename('fullpath')); parentFolder = fileparts(currFolder);
    addpath(genpath(fullfile(currFolder, 'dep')));
    addpath(genpath(fullfile(currFolder, 'plt')));
    addpath(parentFolder);
    params = foragingParams; 
    if isempty(allSubj)
        allSubj = foragingSearchSubj(searchFolder, params);
    end
    allTrlProps = struct([]);
    for subj = allSubj
       subjTrlProps = foragingSubjAnalysis(subj, [], searchFolder, 1, 0);
       allTrlProps = [allTrlProps subjTrlProps]; %#ok<AGROW> 
    end

    clear currFolder parentFolder searchFolder subjTrlProps

    % A ideia básica: 
    % - Cada trial meu é uma unidade experimental para comparar as condições
    %   (sacádica e não sacádica), pois apresentam independência entre si; o 
    %   agrupamento por sujeito é um nível hierárquico adicional.
    % - Temos um design pareado por haver 2 respostas por trial. A diferença
    %   entre respostas num mesmo trial seria a paired difference ou within-
    %   -trial difference. Mas para nós não é relevante.
    % - Por sujeito, conseguimos obter uma subject-level condition difference
    %   que seria a diferença entre d-primes. Isso nos interessa.
    % - Daí nossa estatística seria a média da subject-level condition differ-
    %   ence
    % - Por que preciso preservar relação entre trials e sujeitos? Pois
    %   apesar de a estrutura dos trials ser independente, as respostas
    %   dependem do sujeito, evidentemente.
    % - Se a pergunta for 'há um efeito de condição na população de
    %   sujeitos?', os sujeitos são as unidades independentes para a
    %   inferência a nível populacional
    % - Para cada permutação, para cada sujeito, para cada trial, sortear
    %   ao acaso o rótulo, recalcular a resposta intra-sujeito para cada
    %   condição, a diferença intra-sujeito entre condições e por fim a
    %   estatística do grupo. Com isso, obtenho a distribuição da estatística 
    %   escolhida sob a hipótese nula de troca de rótulos (label
    %   exchangeability) dentro de cada trial
    % - Segundo o chat, eu poderia colocar minha H_0 como: 
    %   Under \(H_0\), conditional on the subject and the experimental design, 
    %   the observed condition assignment is exchangeable according to the 
    %   randomization mechanism used in the experiment.

    %% 1. Efeito pré-sacádico principal
    % ----
    %    d-primes reais por sujeito -> 
    % -> diferenças de d-primes por sujeito ->
    % -> média da diferença
    % ----
    nSubj = numel(allSubj);
    dPerCond = zeros(nSubj, 2);

    for i = 1:nSubj
        subj = allSubj(i);
        subjTrlProps = allTrlProps([allTrlProps.subjNum] == subj);
        [main, ~] = getPSAeffect(subjTrlProps, 0);
        dPerCond(i,:) = [main.nSacc.d, main.sacc.d];
    end
    T = mean(diff(dPerCond, 1, 2));

    % ----
    % Mesmo procedimento que o anterior, mas agora para as permutações
    permT = zeros(1,nPerm);
    for b = 1:nPerm
        permDPerCond = zeros(nSubj, 2);
        permAllTrlProps = foragingPSApermute(allTrlProps);
        for i = 1:nSubj
            subj = allSubj(i);
            subjTrlProps = permAllTrlProps([permAllTrlProps.subjNum] == subj);
            [main, ~] = getPSAeffect(subjTrlProps, 0);
            permDPerCond(i,:) = [main.nSacc.d, main.sacc.d];
        end
        permT(b) = mean(diff(permDPerCond, 1, 2));
    end

    pPSA = plotPSAPermTest(T, permT, params);

    %% 1.a) PSA restrito a distâncias específicas
    % Usar algo como 
%             [closeTrl, farTrl] =plotPSAStimTriangleProps1(probes.pre.pos, probes.sacc.pos, probes.nSacc.pos, mat.drP, 3,0);
%             plotPSAmain(st(closeTrl), mat.drP);
%             plotPSAmain(st(farTrl), mat.drP);
    % E, dependendo, usar as máscaras geradas nas demais análises

    %% 2. Efeito da categoria dos estímulos no desempenho PSA
    % Nesse caso, a hipótese nula seria de que os diferentes níveis dos
    % fatores não interferem no desepenho. Seriam 3 fatores (2x2x2): condição 
    % probe sacc ou não-sacc, categoria do pré-probe e categoria do probe. 
    % Então divide por 4 o número de trials para cada sujeito
    plotPSAcat([allTrlProps.preProbeCat], [allTrlProps.probeCat], [allTrlProps.probeHit], [allTrlProps.nSaccProbeHit], mat.drP)

    %% 3. Efeito do desempenho forrageamento no desempenho PSA
    % ANOVA 2x2, divide em 2 categorias os trials existentes, bem
    % desproporcionais (talvez 5:1). Lembrando que aqui talvez tenha um
    % efeito de confounding com o seguinte, já que a chance de acertos é
    % função da quantidade de vistos

    %% 4. Efeito do tamanho do histórico no desempenho PSA
    % Como pode ir de 2 a 6, ambos inclusos, são 5 categorias diferentes, e
    % emprincípio com a mesma quantidade de trials

    %% 5.a) Efeito da duração do ruído rosa no desempenho PSA
    % Talvez dividir em sacada iniciada durante (i.e., menos de 83 ms) e
    % iniciada depois. ANOVA 2x2

    %% 5.b) Efeito do intervalo sacádico no desempenho PSA
    % Talvez dividir em quartis por sujeito, ou em durações fixas. Também 2
    % fatores

    %% 6. Forrageamento: duração da fixação e acerto

    %% 7. Forrageamento: duração da fixação em função do tamanho do histórico
    % Acredito que tudo funcione melhor se 6 e 7 forem avaliados
    % conjuntamente, o que resulta em 5*(# de níveis de durações), a não
    % ser que use a duração como covariável

    %% 8. Forrageamento: 


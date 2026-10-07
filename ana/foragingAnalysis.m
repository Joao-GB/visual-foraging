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
    allTrlProbes = struct([]);
    for subj = allSubj
       [subjTrlProps, mat, subjTrlProbes] = foragingSubjAnalysis(subj, [], searchFolder, 1, 0);
       allTrlProps = [allTrlProps subjTrlProps]; %#ok<AGROW> 
       allTrlProbes = [allTrlProbes subjTrlProbes]; %#ok<AGROW> 
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

    nSubj = numel(allSubj);
    
    % 1. Efeito pré-sacádico principal
    [pPSA, TPSA, permTPSA] = permAnalysisPSA(allTrlProps, allSubj, nPerm, params);
    
    % 1.a) PSA restrito a distâncias específicas
    %% Veja que talvez faça mais sentido analisar efeito da excentricidade apenas sobre o não sacádico
    modeDist = 1; % 1= direção; 2 = excentricidade; 3 = spotlight          %     [pPSAdist, TPSAdist, permTPSAdist] = permAnalysisPSAdist(allTrlProps, allTrlProbes, allSubj, nPerm, mat, modeDist, params);

    factor2dist = cell(1, nSubj); nLevelDist = 2;
    for subj=1:nSubj
        [closeTrl, farTrl] = plotPSAStimTriangleProps1(allTrlProbes(subj).pre.pos, allTrlProbes(subj).sacc.pos, allTrlProbes(subj).nSacc.pos, mat.drP, modeDist,0);
        factor2dist{subj} = [closeTrl'; farTrl'];
    end

    [dPSAdist, pPSAdist, TPSAdist, permTPSAdist]  = perANOVA_PSA_2byK(allTrlProps, factor2dist, allSubj, nPerm, nLevelDist);
            % Usar algo como 
        %             [closeTrl, farTrl] = plotPSAStimTriangleProps1(allTrlProbes(1).pre.pos, allTrlProbes(1).sacc.pos, allTrlProbes(1).nSacc.pos, mat.drP, 1,0);
        %             plotPSAmain(allTrlProps(closeTrl), mat.drP);
        %             plotPSAmain(allTrlProps(farTrl), mat.drP);
            % E, dependendo, usar as máscaras geradas nas demais análises

    %% 2. Efeito da categoria dos estímulos no desempenho PSA
    %% Aqui só consigo comparar as acurácias, e não os d-primes
    % Nesse caso, a hipótese nula seria de que os diferentes níveis dos
    % fatores não interferem no desepenho. Seriam 3 fatores (2x2x2): condição 
    % probe sacc ou não-sacc, categoria do pré-probe e categoria do probe. 
    % Então divide por 4 o número de trials para cada sujeito
    catMask = calcPSAcat([allTrlProps.preProbeCat], [allTrlProps.probeCat], [allTrlProps.probeHit], [allTrlProps.nSaccProbeHit]); 
%     squeeze([catMask(1,1,:)+catMask(1,2,:);catMask(2,1,:)+catMask(2,2,:)]);
%     catMask = [allTrlProps.preProbeCat];
%     catMask = [catMask; ~catMask];
    nCondCatPreProbe = 2; nCondCatProbe = 2;
    nCondSacc = 2;
    accPerCond = zeros(nSubj,nCondCatPreProbe,nCondCatProbe,nCondSacc);

    for i = 1:nSubj
        subj = allSubj(i);
        for r=1:nCondCatPreProbe
            for c=1:nCondCatProbe
            subjTrlProps = allTrlProps([allTrlProps.subjNum]' == subj & squeeze(catMask(r,c,:)));
            [main, ~] = getPSAeffect(subjTrlProps, 0);
            accPerCond(i,r,c,:) = [main.sacc.table(c,c)/sum(main.sacc.table(c,:)), main.nSacc.table(c,c)/sum(main.nSacc.table(c,:))];
            end
        end
    end
    % Produz as estatísticas de teste
    % (a) Efeito principal condição sacc vs não-sacc: tabela com 2 colunas,
    %     sacc e não sacc, e em cada coluna todas as combinações possíveis
    %     dos demais fatores, CondCatPreProbe e CondCatProbe, em ordem. Na
    %     terceira dimensão, uma 'chapa' por sujeito
    accPerSaccCond = reshape(permute(accPerCond, [3 2 4 1]), ...
            nCondCatProbe*nCondCatPreProbe, nCondSacc, nSubj);
    % Somo ao longo das linhas e das 'chapas', mantendo separados apenas os
    % níveis do fator atualmente em análise
    nEntries = nCondCatProbe * nCondCatPreProbe * nSubj;
    TcatMainSac = sum(squeeze(sum(accPerSaccCond, [1 3])).^2 / nEntries);

    % ----
    permTcatMainSac = zeros(1,nPerm);
    for b = 1:nPerm
        permAccPerCond = zeros(nSubj,nCondCatPreProbe,nCondCatProbe,nCondSacc);
        permAllTrlProps = foragingPSApermute(allTrlProps);
        for i = 1:nSubj
            subj = allSubj(i);
            for r=1:nCondCatPreProbe
                for c=1:nCondCatProbe
                subjTrlProps = allTrlProps([allTrlProps.subjNum]' == subj & squeeze(catMask(r,c,:)));
                [main, ~] = getPSAeffect(subjTrlProps, 0);
                accPerCond(i,r,c,:) = [main.sacc.table(c,c)/sum(main.sacc.table(c,:)), main.nSacc.table(c,c)/sum(main.nSacc.table(c,:))];
                end
            end
        end
        permTcatMainSac(b) = sum(squeeze(sum(permAccPerCond, [1 3])).^2 / nEntries);
    end

    pPSA = plotPSAPermTest(Tmain, permTmain, params);


    % (b) Efeito principal categoria do estímulo
    sum(accPerSaccCond)
    Tmain = sum(sum(dPerCond, 1).^2);

    %% 3. Efeito do desempenho forrageamento no desempenho PSA
    % ANOVA 2x2, divide em 2 categorias os trials existentes, bem
    % desproporcionais (talvez 5:1). Lembrando que aqui talvez tenha um
    % efeito de confounding com o seguinte, já que a chance de acertos é
    % função da quantidade de vistos
    factor2forPerf = cell(1, nSubj); nLevelForPerf = 2;
    for i=1:nSubj
        subj = allSubj(i);
        subjTrlProps = allTrlProps([allTrlProps.subjNum] == subj);
        forProbeHitTrl = [subjTrlProps.forProbeHit] == 1;
        factor2forPerf{i} = [forProbeHitTrl; ~forProbeHitTrl];
    end
    [dPSAforPerf, pPSAforPerf, TPSAforPerf, permTPSAforPerf]  = perANOVA_PSA_2byK(allTrlProps, factor2forPerf, allSubj, nPerm, nLevelForPerf);

    %% 4. Efeito do tamanho do histórico no desempenho PSA
    % Como pode ir de 2 a 6, ambos inclusos, são 5 categorias diferentes, e
    % em princípio com a mesma quantidade de trials
    factor2forHistLen = cell(1, nSubj); nLevelForHistLen = 5;
    for i=1:nSubj
        subj = allSubj(i);
        subjTrlProps = allTrlProps([allTrlProps.subjNum] == subj);
        factor2forHistLen{i} = false(nLevelForHistLen, numel(subjTrlProps));
        for j=1:nLevelForHistLen
            factor2forHistLen{i}(j,:) = [subjTrlProps.forHistLen] == j+1;
        end
    end
    [dPSAforHistLen, pPSAforHistLen, TPSAforHistLen, permTPSAforHistLen]  = perANOVA_PSA_2byK(allTrlProps, factor2forHistLen, allSubj, nPerm, nLevelForHistLen);

    %% 5.a) Efeito da duração do ruído rosa no desempenho PSA
    % Talvez dividir em sacada iniciada durante (i.e., menos de 83 ms) e
    % iniciada depois. ANOVA 2x2
    factor2pinkDur = cell(1, nSubj);
    for i=1:nSubj
        subj = allSubj(i);
        subjTrlProps = allTrlProps([allTrlProps.subjNum] == subj);
        pinkDurMask = plotPSApinkDur(subjTrlProps, mat, [], []);
        factor2pinkDur{i} = pinkDurMask;
    end
    nLevelPinkDur = size(pinkDurMask,1);
    [dPSApinkDur, pPSApinkDur, TPSApinkDur, permTPSApinkDur]  = perANOVA_PSA_2byK(allTrlProps, factor2pinkDur, allSubj, nPerm, nLevelPinkDur);



    %% 5.b) Efeito do intervalo sacádico no desempenho PSA
    % Talvez dividir em quartis por sujeito, ou em durações fixas. Também 2
    % fatores

    %% 6. Forrageamento: duração da fixação e acerto

    %% 7. Forrageamento: duração da fixação em função do tamanho do histórico
    % Acredito que tudo funcione melhor se 6 e 7 forem avaliados
    % conjuntamente, o que resulta em 5*(# de níveis de durações), a não
    % ser que use a duração como covariável

    %% 8. Forrageamento: 


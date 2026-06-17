function [trl, eyeData, eventTimes] = foragingTrlProps(mat, edf)
    currFolder = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(currFolder, 'dep')));

    N = mat.tkP.nBlocks*mat.tkP.nTrials;
    trl = getTrlStruct(N);

    % Mensagens do arquivo
    messages = {edf.FEVENT(:).message}; messages(cellfun(@isempty, messages)) = {''};

    events = {edf.FEVENT(:).codestring};
    eventTimes  = [edf.FEVENT(:).sttime; edf.FEVENT(:).entime];

    % Extrai o que importa quanto à trajetória do olho
    Eye = mat.tkP.Eye;
    eyeData.coord = [edf.FSAMPLE.gx(Eye, :);edf.FSAMPLE.gy(Eye, :)];
    eyeData.time = edf.FSAMPLE.time;
    fixEvents = {'STARTFIX', 'ENDFIX'};
    
    %% Encontra o início das sessões experimentais
    sesOnMsg = sprintf(mat.prm.msg.on.ses{1}, mat.prm.msg.suffix{3});
    expSesOnIdx = [find(contains(messages, sesOnMsg))];
    
    if ~isempty(expSesOnIdx)
        auxIdx = [expSesOnIdx length(messages)];
        expSesOffIdx = zeros(1, length(expSesOnIdx));
        for i=1:length(expSesOnIdx)
            offAuxIdx = [find(contains(messages(auxIdx(i):auxIdx(i+1)), mat.prm.msg.off.ses{3}))] + expSesOnIdx(i) - 1;
            if length(offAuxIdx) ~= 1
                error('O arquivo deve conter exatamente uma sessão experimental concluída!')
            else
                expSesOffIdx(i) = offAuxIdx;
            end
        end
        expSesOnIdx(expSesOffIdx == 0)  = [];
        expSesOffIdx(expSesOffIdx == 0) = [];

        %% Encontra a última sessão experimental completa
        sesCpMsg = sprintf(mat.prm.msg.off.ses{1}, mat.prm.msg.suffix{3});
        auxIdx = [find(strcmp(messages(expSesOffIdx), sesCpMsg), 1, "last")];
        if isempty(auxIdx)
            warning('Não há sessões experimentais concluídas, usando a última interrompida');
            sesCpMsg = sprintf(mat.prm.msg.off.ses{2}, mat.prm.msg.suffix{3});
            auxIdx = min([length(messages), [find(strcmp(messages(expSesOffIdx), sesCpMsg), 1, "last")]]);
        end
    
        sesOnIdx  = expSesOnIdx(auxIdx);
        sesOffIdx = expSesOffIdx(auxIdx);
        sesLimIdx = [sesOnIdx; sesOffIdx];
        clear sesOnMsg offAuxIdx expSesOnIdx expSesOffIdx sesOnIdx sesOffIdx auxIdx i
    
        %% Encontra o onset e offset dos trials
            % Assumo que meu programa retorna fim para todos os trials, inclusive
            % se não for bom
        trlOnIdx = find(contains(messages(sesLimIdx(1):sesLimIdx(2)), mat.prm.msg.on.trl{2})) + sesLimIdx(1) - 1;
        trlOffIdx = zeros(1, length(trlOnIdx)); auxIdx = [trlOnIdx sesLimIdx(2)];
        for i=1:length(trlOnIdx)
            aux = find(contains(messages(auxIdx(i):auxIdx(i+1)), mat.prm.msg.off.trl{2}), 1) + trlOnIdx(i) - 1;
            if isempty(aux)  %% Apenas para qualquer erro que aconteça de eu não pausar direito 
                trlOffIdx(i) = trlOnIdx(i+1);
            else
                trlOffIdx(i) = aux;
            end
        end
        trlLimIdx = [trlOnIdx; trlOffIdx];
        trlLimAbs = reshape(eventTimes(1, trlLimIdx), size(trlLimIdx));
        trlDur = double(trlLimAbs(2,:)-trlLimAbs(1,:))/1000;

        nTrl = size(trlLimIdx, 2);
    
        %% Encontra os trials descartados por pausa, timeout etc.
        badTrl = zeros(1, nTrl);
        for i=1:nTrl
            currTrlAllIdx = trlLimIdx(1,i):trlLimIdx(2,i);
            badCond = strcmp(messages(currTrlAllIdx), mat.prm.msg.err.trl{1}) | ... % Cobre efeitos de pausa
                ... strcmp(messages(currTrlIdx), mat.prm.msg.err.trl{2}) | ...   % Não incluo efeito de abort idx para manter nTrl 
                ...                                                              % por bloco. Essa info já está em trialOrder(2,:,:)
                strcmp(messages(currTrlAllIdx), mat.prm.msg.err.P1) | ...
                strcmp(messages(currTrlAllIdx), mat.prm.msg.err.P2);
            if any(badCond), badTrl(i) = true; end
        end
    
        %% Encontra os trials que estão em blocos descartados e os adiciona aos descartados
        sesItMsg = sprintf(mat.prm.msg.off.ses{2}, mat.prm.msg.suffix{3});
        blkOnIdx = find(contains(messages(sesLimIdx(1):sesLimIdx(2)), mat.prm.msg.on.blk{2})) + sesLimIdx(1) - 1;
        blkOffIdx = zeros(1, length(blkOnIdx));
        nBlk = length(blkOnIdx);
        auxIdx = [blkOnIdx sesLimIdx(2)]; badBlk = zeros(1, nBlk);
        for i=1:nBlk
            currBlkIdx = auxIdx(i):auxIdx(i+1);
            blkOffIdx(i) = find(contains(messages(currBlkIdx), mat.prm.msg.off.blk{2})) + auxIdx(i) - 1;
            badCond = strcmp(messages(currBlkIdx), mat.prm.msg.err.blk) | ...
                strcmp(messages(currBlkIdx), sesItMsg);
            if any(badCond), badBlk(i) = true; end
        end
        blkIdx = [blkOnIdx; blkOffIdx];
        badblkIdx = find(badBlk);
        for i=1:length(badblkIdx)
            badTrl(trlLimIdx(1,:) >= auxIdx(badblkIdx(i)) & trlLimIdx(1,:) <= auxIdx(badblkIdx(i)+1)) = true;
        end
        blkIdx(:, badblkIdx) = [];
    
        %% Deleta os trials descartados
        trlLimIdx(:, logical(badTrl)) = [];
        nTrl = size(trlLimIdx, 2);
        if nTrl ~=mat.tkP.nTrials*mat.tkP.nBlocks
            warning('O número de trials não é compatível com o número esperado')
        end
    
        clear aux sesItMsg sesCpMsg auxIdx badblkIdx badBlk currBlkIdx blkOnIdx...
                goodTrl badTrl badCond currTrlAllIdx trlOnIdx trlOffIdx nBlk i
        % Lembre que a fila de estímulos vistos é indexada na ordem em que
        % os trials foram apresentados, então stimsQueue{1} é o 1o trial do 
        % 1o bloco, stimsQueue{mat.tkP.nTrials+1) é o 1o do 2o bloco etc.,
        % por mais que, e.g., o índice do 1o trial tenha sido 18.
        stimsQueue = reshape(mat.results.seenStimsQueue', 1, []);

        %% Extrai todas as fixações segundo as mensagens que escrevi no Eyelink
        % O limIdx é da forma: 
        % [início da fix; 
        % fim da fix;
        % duração da fix;         (exata, com base na 
        % % índice do estímulo visto segundo as matrizes de results;
        % visita a novo estímulo; (1 = sim)
        % fases em que ocorre;    (2 = P2, 3 = P3, 4 = PM, 5 = P2+P3, 7 = P3+P4)
        % término a força]        (1 = sim)

        b = 1;

        for i=1:nTrl
            t = mod(i-1, mat.tkP.nTrials)+1;
            % Um trial não é mantido apenas se anteriormente não era bom (o
            % que acontece quando aborta o índice por excesso de repetições),
            % ou porque alguma das fases (p1 a pM) não apresenta a quantidade
            % necessária de estímulos
            trl(i).trlKeep = mat.results.trialOrder(2,t,b);
            if trl(i).trlKeep
                %% Gerais
                trl(i).trlLimIdx = trlLimIdx(:,i);
                trl(i).trlLimAbs = trlLimAbs(:,i);
                trl(i).trlDur = trlDur(i);
                trl(i).trlIdx = mat.results.trialOrder(1,t,b);
                trl(i).blkIdx = b;
                trl(i).tgtOri = mat.results.targetOri(b);
                trl(i).nTgt   = mat.results.nTs(b,trl(i).trlIdx);
                trl(i).tgtIdx = find(mat.results.orientation(:,trl(i).trlIdx, b) == trl(i).tgtOri);
                trl(i).tgtPos = mat.results.stimCenters(:,trl(i).tgtIdx,trl(i).trlIdx,b);
                trl(i).stmOri = mat.results.orientation(:,trl(i).trlIdx, b);
                trl(i).stmPos = mat.results.stimCenters(:,:,trl(i).trlIdx,b);

                %% Fase 1
                trl(i).p1Pos  = mat.results.fixCenters(:, trl(i).trlIdx, b);

                % Todos os índices do trial
                currTrlAllIdx = trlLimIdx(1,i):trlLimIdx(2,i);

                [pLI, pLA, pD] = getPhaseLims(messages, eventTimes(1,:), currTrlAllIdx, mat.prm.msg.on.P1, mat.prm.msg.on.P2);
                trl(i).p1LimIdx = pLI;
                trl(i).p1LimAbs = pLA;
                trl(i).p1Dur    = pD;

                [fLI, fLA, fTI, fD, fEP, fCV, detail] = getFixLims(events, eventTimes, eyeData.time, trl(i).p1LimIdx(2), fixEvents, trl(i).p1LimIdx, eyeData.coord);
                if isempty(fLI), trl(i).trlKeep = 0; 
                    continue; 
                end
                trl(i).p1FixLimIdx = fLI;
                trl(i).p1FixLimAbs = fLA;
                trl(i).p1FixTimeIdx = fTI;
                trl(i).p1FixDur    = fD;
                trl(i).p1FixPos    = fEP;
                trl(i).p1FixSpread = fCV;
                trl(i).p1FixDetail = detail;
                
               
                %% Fase 2
                [pLI, pLA, pD] = getPhaseLims(messages, eventTimes(1,:), currTrlAllIdx, mat.prm.msg.on.P2, mat.prm.msg.on.P3);
                trl(i).p2LimIdx = pLI;
                trl(i).p2LimAbs = pLA;
                trl(i).p2Dur    = pD;

                % Encontra os estímulos dentro do trial
                p2StmOnIdx = find(strcmp(messages(trl(i).p2LimIdx(1):trl(i).p2LimIdx(2)), mat.prm.msg.on.stm{1})|strcmp(messages(trl(i).p2LimIdx(1):trl(i).p2LimIdx(2)), mat.prm.msg.on.stm{2})) + trl(i).p2LimIdx(1) - 1;
                auxIdx = [p2StmOnIdx trlLimIdx(2,i)];
                p2StmOffIdx = arrayfun(@(j) find(strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.off.stm{1})|strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.off.stm{3})|strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.off.stm{2}), 1) + auxIdx(j) - 1,  1:length(p2StmOnIdx));
                aux = strcmp(messages(p2StmOffIdx), mat.prm.msg.off.stm{2});
                p2StmOnIdx(aux) = []; p2StmOffIdx(aux) = [];
                p2StmIdx = [p2StmOnIdx; p2StmOffIdx];
%                 
                [fLI, fLA, fTI, fD, fEP, fCV, detail] = getFixLims(events, eventTimes, eyeData.time, trl(i).p2LimIdx, fixEvents, p2StmIdx, eyeData.coord);
                if isempty(fLI), trl(i).trlKeep = 0;
                    continue;
                end
                trl(i).p2FixLimIdx = fLI;
                trl(i).p2FixLimAbs = fLA;
                trl(i).p2FixTimeIdx = fTI;
                trl(i).p2FixDur    = fD;
                trl(i).p2FixPos    = fEP;
                trl(i).p2FixSpread = fCV;
                trl(i).p2FixDetail = detail;

                trl(i).p2NFix = numel(p2StmOnIdx);
                trl(i).p2FixIdx = stimsQueue{i}(1,1:trl(i).p2NFix);
                trl(i).p2NFixNoRep = numel(unique(trl(i).p2FixIdx));
                trl(i).p2OnlyNFix = trl(i).p2NFix - 1;
                trl(i).p2OnlyNFixNoRep = trl(i).p2NFixNoRep - 1;
                trl(i).p2OnlyFixIdx = trl(i).p2FixIdx(1:end-1);
                trl(i).p2FixOnTgt = ismember(trl(i).p2FixIdx, trl(i).tgtIdx);
                
                


                if trl(i).p2NFixNoRep ~= mat.results.modTimes(b, trl(i).trlIdx)
                    disp('Há algo de ERRADO com fase 2!'); 
                end

                trl(i).p2Pos = trl(i).stmPos(:,trl(i).p2FixIdx);

                trl(i).p2PosDist = pointsDist([trl(i).p1Pos trl(i).p2Pos]);
                trl(i).p2FixDist = pointsDist([trl(i).p1FixPos trl(i).p2FixPos]);


                %% Fase 3 (ruído rosa)
                [pLI, pLA, pD] = getPhaseLims(messages, eventTimes(1,:), currTrlAllIdx, mat.prm.msg.on.P3, mat.prm.msg.off.P3);
                trl(i).p3LimIdx = pLI;
                trl(i).p3LimAbs = pLA;
                trl(i).p3Dur    = pD;

                % Encontra os estímulos dentro da fase 3 do trial
                    % Com certeza vai ter um Off, seja por movimento
                    % ocular ou por fim da fase 3
                p3StmOffIdx = find(strcmp(messages(trl(i).p3LimIdx(1):trl(i).p3LimIdx(2)), mat.prm.msg.off.stm{1})|strcmp(messages(trl(i).p3LimIdx(1):trl(i).p3LimIdx(2)), mat.prm.msg.off.stm{3})) + trl(i).p3LimIdx(1) - 1;
                auxIdx = [trlLimIdx(1,i) p3StmOffIdx];
                p3StmOnIdx = arrayfun(@(j) find(strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.on.stm{1})|strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.on.stm{2}), 1,'last') + auxIdx(j) - 1,  1:length(p3StmOffIdx));
                p3StmIdx = [p3StmOnIdx; p3StmOffIdx];

                % MAS não necessariamente haverá início ou fim de uma fixação, o mais provável, se
                % a pessoa executar a tarefa conforme o previsto, é que haja ao menos o fim de uma fixação
                % (nunca haverá apenas um início, pois por construção a fase 3 inicia quando a última fixação
                % se inicia). Se não houver, pode ser que a sacada foi feita um pouco depois
                if p3StmIdx(2,end) <= p2StmIdx(2,end)
                    [fLI, fLA, fTI, fD, fEP, detail] = getFixLims(events, eventTimes, eyeData.time, [trl(i).p2FixLimIdx(1,end); trl(i).p2FixLimIdx(2,end)], fixEvents, [p3StmIdx(1); trl(i).p2FixLimIdx(2,end)], eyeData.coord, 0);
                else
                    [fLI, fLA, fTI, fD, fEP, fCV, detail] = getFixLims(events, eventTimes, eyeData.time, trl(i).p3LimIdx, fixEvents, p3StmIdx, eyeData.coord, 0);
                end
                if isempty(fLI)
                    trl(i).trlKeep = 0; 
                    continue;
                end
                trl(i).p3FixLimIdx = fLI;
                trl(i).p3FixLimAbs = fLA;
                trl(i).p3FixTimeIdx = fTI;
                trl(i).p3FixDur    = fD;
                trl(i).p3FixPos    = fEP;
                trl(i).p3FixSpread = fCV;
                trl(i).p3FixDetail = detail;

                trl(i).p3NFix = numel(p3StmOnIdx);
                trl(i).p3OnlyNFix = trl(i).p3NFix - 1;
                trl(i).p3FixIdx = stimsQueue{i}(1,(trl(i).p2OnlyNFix+1):(trl(i).p2OnlyNFix+1+trl(i).p3OnlyNFix));
                trl(i).p3OnlyFixIdx = setdiff(trl(i).p3FixIdx, trl(i).p2FixIdx);
                trl(i).p3OnlyFixOnTgt = ismember(trl(i).p3OnlyFixIdx, trl(i).tgtIdx);
                trl(i).p3p2FixIdx = intersect(trl(i).p3FixIdx, trl(i).p2FixIdx);

                % O período durante o qual foi visto o ruído rosa é o
                % mínimo entre a duração da fase 3 e a diferença entre seu
                % início e o fim da última fixação durante essa fase
                trl(i).pinkNoiseDur = min((double(trl(i).p3FixLimAbs(2, end)) - double(trl(i).p3LimAbs(1)))/1000, trl(i).p3Dur);

                % O ideal é que sejam valores negativos, pois indica que o
                % movimento ocular que interrompeu o ruído rosa e o ruído 
                % acaba ligeiramente depois 
                trl(i).p3FixOffsetDelay = (double(trl(i).p3FixLimAbs(2, end)) - double(trl(i).p3LimAbs(2)))/1000;
                

                if trl(i).p3OnlyNFix > 0
                    disp('Há algo de ERRADO com fase 3!');
                end

                trl(i).p3Pos = trl(i).stmPos(:,trl(i).p3FixIdx);

                % Se eu quiser verificar se houve mudança da última fixação de
                % P2 não em P3 para P3, vai ser a mesma que entre os dois últimos de P2
%                 auxIdx = numel(setdiff(trl(i).p2FixIdx, trl(i).p3FixIdx));
%                 trl(i).p3PosDist = pointsDist([trl(i).p2Pos(:,auxIdx) trl(i).p3Pos]);
%                 trl(i).p3FixDist = pointsDist([trl(i).p2FixPos(:,auxIdx) trl(i).p3FixPos]);

                trl(i).p3PosDist = pointsDist([trl(i).p2Pos(:,end) trl(i).p3Pos]);
                trl(i).p3FixDist = pointsDist([trl(i).p2FixPos(:,end) trl(i).p3FixPos]);


                %% Pós-modificação
                [pLI, pLA, pD] = getPhaseLims(messages, eventTimes(1,:), currTrlAllIdx, mat.prm.msg.on.PM, mat.prm.msg.off.PM);
                trl(i).pMLimIdx = pLI;
                trl(i).pMLimAbs = pLA;
                trl(i).pMDur    = pD;

                % Encontra os estímulos dentro da PM do trial
                pMStmOffIdx = find(strcmp(messages(trl(i).pMLimIdx(1):trl(i).pMLimIdx(2)), mat.prm.msg.off.stm{1})|strcmp(messages(trl(i).pMLimIdx(1):trl(i).pMLimIdx(2)), mat.prm.msg.off.stm{4})) + trl(i).pMLimIdx(1) - 1;
                auxIdx = [trlLimIdx(1,i) pMStmOffIdx];
                pMStmOnIdx = arrayfun(@(j) find(strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.on.stm{1})|strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.on.stm{2}), 1,'last') + auxIdx(j) - 1,  1:length(pMStmOffIdx));
                pMStmIdx = [pMStmOnIdx; pMStmOffIdx];
                

                % Espero que na PM haja preferencialmente o início de uma
                % fixação (se P3 encerrada por movimento ocular); pode ser
                % que inclua um pouco do fim da anterior, mas se não tiver
                % início é inútil
                [fLI, fLA, fTI, fD, fEP, fCV, detail] = getFixLims(events, eventTimes, eyeData.time, trl(i).pMLimIdx, fixEvents, pMStmIdx, eyeData.coord);
                if isempty(fLI)
                    if ~isempty(pMStmIdx)
                        [fLI, fLA, fTI, fD, fEP, fCV, detail] = getFixLims(events, eventTimes, eyeData.time, [trl(i).p3LimIdx(1); trl(i).pMLimIdx(2)], fixEvents, pMStmIdx, eyeData.coord);
                    end
                    if isempty(fLI)
                        trl(i).trlKeep = 0; 
                        continue;
                    end
                end
                if all(fLI(1,:) < trl(i).pMLimIdx(1))
                    trl(i).trlKeep = 0; 
                    continue;
                end
                trl(i).pMFixLimIdx = fLI;
                trl(i).pMFixLimAbs = fLA;
                trl(i).pMFixTimeIdx = fTI;
                trl(i).pMFixDur    = fD;
                trl(i).pMFixPos    = fEP;
                trl(i).pMFixSpread = fCV;
                trl(i).pMFixDetail = detail;

                % Encontra o primeiro STARTFIX dentro do PM
                auxIdx = find(cellfun(@(x) strcmp(x, fixEvents(1)), events(trl(i).pMLimIdx(1):trl(i).pMLimIdx(2))), 1) + trl(i).pMLimIdx(1) - 1;
                auxIdx = find(cellfun(@(x) strcmp(x, fixEvents(2)), events(auxIdx:trl(i).trlLimIdx(2))), 1) + auxIdx - 1;
                if ~isempty(auxIdx)
                    trl(i).pMFirstFixPos = mean(eyeData.coord(:,find(eyeData.time == eventTimes(1,auxIdx), 1):find(eyeData.time == eventTimes(2,auxIdx), 1)),2);
                end
                trl(i).pMNFix = numel(pMStmOnIdx);
                trl(i).pMp3FixIdx = find(ismember(trl(i).pMFixLimIdx(1,:), trl(i).p3FixLimIdx(1,:)));
                trl(i).pMOnlyNFix = trl(i).pMNFix - numel(trl(i).pMp3FixIdx);
                trl(i).pMFixIdx = stimsQueue{i}(1,(numel(trl(i).p3p2FixIdx)+trl(i).p2OnlyNFix+trl(i).p3OnlyNFix+1):(numel(trl(i).p3p2FixIdx)+trl(i).p2OnlyNFix+trl(i).p3OnlyNFix+trl(i).pMOnlyNFix));
                trl(i).pMOnlyFixIdx = setdiff(trl(i).pMFixIdx, [trl(i).p3FixIdx]);
                trl(i).pMOnlyFixOnTgt = ismember(trl(i).pMOnlyFixIdx, trl(i).tgtIdx);
                if isempty(trl(i).pMOnlyFixIdx)
                    trl(i).trlKeep = 0; 
                    continue;
                end


                if trl(i).pMOnlyNFix == 0
                    disp('Há algo de ERRADO com PM!'); 
                end


                trl(i).pMPos = trl(i).stmPos(:,trl(i).pMFixIdx);


                trl(i).pMPosDist = pointsDist([trl(i).p3Pos(:,end) trl(i).pMPos]);
                trl(i).pMFixDist = pointsDist([trl(i).p3FixPos(:,end) trl(i).pMFixPos]);

                %% Fase 4

                [pLI, pLA, pD] = getPhaseLims(messages, eventTimes(1,:), currTrlAllIdx, mat.prm.msg.on.P4, mat.prm.msg.off.P4);
                trl(i).p4LimIdx = pLI;
                trl(i).p4LimAbs = pLA;
                trl(i).p4Dur    = pD;

                trl(i).p4IdxAsk = mat.results.trialFeedback{b,t}(1,:);
                trl(i).p4IdxStatus = mat.results.trialFeedback{b,t}(2,:);
                trl(i).p4IdxFB = mat.results.trialFeedback{b,t}(3,:);

                trl(i).p4PosDist = pointsDist(trl(i).pMPos, trl(i).stmPos(:,trl(i).p4IdxAsk));
                trl(i).p4PosDistm1 = trl(i).p4PosDist(:, trl(i).p4IdxStatus == -1);
                trl(i).p4PosDist0 = trl(i).p4PosDist(:, trl(i).p4IdxStatus == 0);
                trl(i).p4PosDist1 = trl(i).p4PosDist(:, trl(i).p4IdxStatus == 1);

                trl(i).p3ToP4PosDist = pointsDist(trl(i).stmPos(:,trl(i).p3p2FixIdx), trl(i).stmPos(:,trl(i).p4IdxAsk));
                trl(i).p3ToP4PosDistm1 = trl(i).p3ToP4PosDist(:, trl(i).p4IdxStatus == -1);
                trl(i).p3ToP4PosDist0  = trl(i).p3ToP4PosDist(:, trl(i).p4IdxStatus == 0);
                trl(i).p3ToP4PosDist1  = trl(i).p3ToP4PosDist(:, trl(i).p4IdxStatus == 1);


%                 trl(i).p4AnsDur = 1;

                trl(i).p4AnsExp = ismember(trl(i).p4IdxAsk, trl(i).tgtIdx);
                trl(i).p4AnsExpm1 = trl(i).p4AnsExp(trl(i).p4IdxStatus == -1);
                trl(i).p4AnsExp0 = trl(i).p4AnsExp(trl(i).p4IdxStatus == 0);
                trl(i).p4AnsExp1 = trl(i).p4AnsExp(trl(i).p4IdxStatus == 1);

                trl(i).p4AnsObs = xor(trl(i).p4AnsExp, ~trl(i).p4IdxFB);
                trl(i).p4AnsObsm1 = trl(i).p4AnsObs(trl(i).p4IdxStatus == -1);
                trl(i).p4AnsObs0 = trl(i).p4AnsObs(trl(i).p4IdxStatus == 0);
                trl(i).p4AnsObs1 = trl(i).p4AnsObs(trl(i).p4IdxStatus == 1);

                trl(i).p4AnsFBm1 = trl(i).p4IdxFB(trl(i).p4IdxStatus == -1);
                trl(i).p4AnsFB0  = trl(i).p4IdxFB(trl(i).p4IdxStatus == 0);
                trl(i).p4AnsFB1  = trl(i).p4IdxFB(trl(i).p4IdxStatus == 1);

                trl(i).p4StatusHas0and1 = all(ismember([0 1], trl(i).p4IdxStatus));
                trl(i).p4StatusHas0or1  = any(ismember([0 1], trl(i).p4IdxStatus));
                
            end
            if t == mat.tkP.nTrials, b = b+1; end
        end

    

    else
        disp('ERRO: não há sessões experimentais no arquivo!')
    end
        
end

function [pLI, pLA, pD] = getPhaseLims(m, st, tAI, mOn, mOff)
    pLI = [find(strcmp(m(tAI), mOn),1); find(strcmp(m(tAI), mOff),1)] + tAI(1) - 1;
    pLA = reshape(st(pLI), [2,1]);
    pD = double(pLA(2)-pLA(1))/1000;
end

function [fLI, fLA, fTI, fD, fEP, fCV, detail] = getFixLims(e, et, times, pLims, fEv, fLims, eyePos, mode)
    % O modo indica se a fixação será procurada com relação ao fim ou ao
    % início do intervalo
    fLI=[]; fLA=[]; fTI=[]; fD=[]; fEP=[]; fCV = []; detail=[];
    %% 1 se quiser lidar com início de fixação; 0 para fins (P3, pois fixação começou na P2)
    if nargin < 8, mode = 1; end
    if mode
        if numel(pLims) == 1 % Procura a fixação mais próxima a começar antes
            aux = find(strcmp(e(1:pLims), fEv{1}), 1, "last");
        else                % Procura todas as fixações no intervalo, 
                aux = find(strcmp(e(pLims(1):pLims(2)), fEv{1})) + pLims(1) - 1;
        end
        N = numel(aux);
        if N == 0, return;
        else
            aux2 = find(strcmp(e(aux(1):end), fEv{2}),N) + aux(1) - 1;
            if numel(aux) ~= numel(aux2), return;
            else
                fLI = [aux; aux2];
            end
        end
    else
        if numel(pLims) == 1 % Procura a fixação mais próxima a terminar antes
            aux = find(strcmp(e(1:pLims), fEv{2}), 1, "last");
        else                % Procura todas as fixações no intervalo, 
            aux = find(strcmp(e(pLims(1):pLims(2)), fEv{2})) + pLims(1) - 1;
        end
        N = numel(aux);
        if N == 0, return;
        else
            aux2 = find(strcmp(e(1:aux(N)), fEv{1}),N, "last");
            if numel(aux2)~=numel(aux), return;
            else
                fLI = [aux2; aux];
            end
        end
    end


    fLA = reshape(et(:, fLI(2,:)), [2,N]);
    fTI = [arrayfun(@(x) find(times==x,1), fLA(1,:)); arrayfun(@(x) find(times==x,1), fLA(2,:))];
    fD = double(fLA(2,:)-fLA(1,:))/1000;

    fEP = zeros(2, size(fTI,2)); fCV = zeros(2, 2, size(fTI,2));
    for i=1:size(fTI,2)
        [fEP(:,i), fCV(:,:,i)] = fixStats(eyePos(:,fTI(1,i):fTI(2,i)));
    end

    % Agrega as fixações conforme os limites esperados (e.g., se houver
    % duas fixações dentro de um mesmo estímulo, considera como uma só

    nFix = size(fLims, 2);
    fLims(1,1) = 1; fLims(2,end) = numel(e);
    auxfLI = {}; auxfLA = {}; auxfTI = {}; auxfD = {}; auxfEP = {}; auxfCV = {};
    rfLI = zeros(2,nFix); rfLA = zeros(2,nFix); rfTI = zeros(2,nFix); rfD = zeros(1,nFix); rfEP = zeros(2,nFix); rfCV = zeros(2,2,nFix);

    for i=1:nFix
        % talvez não faça sentido incluir todas que comecem ou terminem
        idx = find((fLI(1,:) >= fLims(1,i) & fLI(1,:) <= fLims(2,i)) | (fLI(2,:) >= fLims(1,i) & fLI(2,:) <= fLims(2,i)));
        if isempty(idx)
            idx = find((fLI(1,:) <= fLims(1,i) & fLI(2,:) >= fLims(1,i)) | (fLI(1,:) <= fLims(2,i) & fLI(2,:) >= fLims(2,i)));
        end

        auxfLI = [auxfLI {fLI(:,idx)}]; %#ok<*AGROW> 
        auxfLA = [auxfLA {fLA(:,idx)}];
        auxfTI = [auxfTI {fTI(:,idx)}];
        auxfD = [auxfD {fD(idx)}];
        auxfEP = [auxfEP {fEP(:,idx)}];
        auxfCV = [auxfCV {fCV(:,:,idx)}];

        rfLI(:,i) = [fLI(1,idx(1));fLI(2,idx(end))];
        rfLA(:,i) = [fLA(1,idx(1));fLA(2,idx(end))];
        rfTI(:,i) = [fTI(1,idx(1));fTI(2,idx(end))];
        rfD(i) = double(rfLA(2,i)-rfLA(1,i))/1000;
        auxD = auxfD{end}./sum(auxfD{end}(:));
        rfEP(:,i) = auxfEP{end}*auxD';

        K = length(idx);
        n = auxfD{end}(:)';
        outfCV = zeros(2,2); outfEP = rfEP(:,i);
        
        for j = 1:K
            
            mu_i = auxfEP{end}(:,j);
            Sigma_i = auxfCV{end}(:,:,j);
            
            diff = mu_i - outfEP;
            
            outfCV = outfCV + ...
                (n(j)-1)*Sigma_i + ...
                n(j)*(diff*diff');
        end
        rfCV(:,:,i) = outfCV;
    end

    detail.fLI = auxfLI; detail.fLA = auxfLA;
    detail.fTI = auxfTI; detail.fD  = auxfD; detail.fEP = auxfEP;
    detail.fCV = auxfCV;
    fLI = rfLI; fLA = rfLA;
    fTI = rfTI; fD  = rfD; fEP = rfEP; fCV = rfCV;

end

function d = pointsDist(p, q)
    if nargin == 1
        dxdy = diff(p, [], 2);

        r = vecnorm(dxdy);
        a = atan2d(dxdy(2,:), dxdy(1,:)); % Ângulos em graus
    % Se tiver 2 argumentos, calcula a distância entre os pontos de p e os
    % pontos de q par a par, resultando em tamanho 2 X (numel(p) x numel(q))
    else
        dx = p(1,:)' - q(1,:); dx = dx(:)';
        dy = p(2,:)' - q(2,:); dy = dy(:)';
        
        r = sqrt(dx.^2 + dy.^2);
        a = atan2d(dy, dx);
    end
    d = [r; a];
end


function [mu, Sigma] = fixStats(X)
% X must be 2 x N

Xn = X';           % N x 2
mu = mean(Xn,1)';  % 2 x 1
Sigma = cov(Xn);   % 2 x 2

end


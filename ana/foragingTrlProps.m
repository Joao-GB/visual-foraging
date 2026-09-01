function [trl, eyeData, eventLimClk] = foragingTrlProps(mat, edf, sesStr, subj, expIdx)
    if nargin < 5, expIdx = 4; end

    if expIdx == 2
        N = numel(mat.prm.allOri)*mat.tkP.t1.tkP.nTrainTrials;
        mat.results = mat.tkP.t1;
        actTrl = mat.tkP.t1.tkP.nTrainTrials;
    else        
        N = mat.tkP.nBlocks*mat.tkP.nTrials;
        actTrl = mat.tkP.nTrials;
    end
    trl    = getTrlStruct(N);

    fprintf('\n--------------\nLEMBRE-SE: tanto o eixo y dos estímulos como da posição ocular serão invertidos\n--------------\n')
    % Para uma versão errada da tarefa, em que isSaccSeen estava transposta
    if size(mat.results.isSaccSeen,1) == size(mat.results.isP3earlyStop,2)
        mat.results.isSaccSeen = mat.results.isSaccSeen(:, 1:mat.tkP.nBlocks)';
    end

    %% Importa mensagens e eventos do arquivo edf
    messages = {edf.FEVENT(:).message}; messages(cellfun(@isempty, messages)) = {''};
%     events = {edf.FEVENT(:).codestring};
    eventLimClk  = [edf.FEVENT(:).sttime; edf.FEVENT(:).entime];

    % Extrai o que importa quanto à trajetória do olho
    Eye = mat.tkP.Eye;
    eyeData.coord = [edf.FSAMPLE.gx(Eye, :);
                    -edf.FSAMPLE.gy(Eye, :) + mat.dpP.winRect(4)];
    eyeData.time = edf.FSAMPLE.time;
    eyeData.pupil = [];
    if ~isfield(mat.prm, 'fs'), mat.prm.fs = 1000; end 
    eyeData.fs    = mat.prm.fs;
    eyeData.sessionID = sesStr;
%     fixEvents = {'STARTFIX', 'ENDFIX'};
    if numel(mat.prm.msg.on.stm) < 3, mat.prm.msg.on.stm{3} = 'STIM ONSET'; end
    stmMsgs = {mat.prm.msg.on.stm{3}, mat.prm.msg.off.stm{1}};
    
    %% Encontra limites da sessão experimental
    sesLimIdx = getSesLims(mat, messages, expIdx);
    if isempty(sesLimIdx)
        error('ERRO: não há sessões experimentais no arquivo!');
    end

    %% Encontra limites dos trials
    [trlLimIdx, trlLimAbs, trlDur] = getTrlLims(mat, messages, sesLimIdx, eventLimClk);
    nTrl = size(trlLimIdx, 2);

    %% Exclui os trials ruins
    [badTrl, nTrl, ~] = delBadTrl(mat, messages, sesLimIdx, trlLimIdx, nTrl);
    trlLimIdx(:,badTrl) = [];
    trlLimAbs(:,badTrl) = [];
    trlDur(badTrl) = [];

    if nTrl ~= N
        warning('O número de trials não é compatível com o número esperado')
    end

    %% Extrai as propriedades de cada trial
    screenDist  = mat.prm.screenDist;
    screenWidth = mat.dpP.monitorW_mm/10;
    screenRes   = mat.dpP.winRect(3:4);
    screenCenter= mat.dpP.winCenter';

    % Lembre que a fila de estímulos vistos é indexada na ordem em que
    % os trials foram apresentados, então stimsQueue{1} é o 1o trial do 
    % 1o bloco, stimsQueue{mat.tkP.nTrials+1) é o 1o do 2o bloco etc.,
    % por mais que, e.g., o índice do 1o trial tenha sido 18.
    stimsQueue = reshape(mat.results.seenStimsQueue', 1, []);

    % O limIdx é da forma: 
    % [início da fix; 
    % fim da fix;
    % duração da fix;         (exata, com base na 
    % % índice do estímulo visto segundo as matrizes de results;
    % visita a novo estímulo; (1 = sim)
    % fases em que ocorre;    (2 = P2, 3 = P3, 4 = PM, 5 = P2+P3, 7 = P3+P4)
    % término a força]        (1 = sim)

    b = 1;

    hasProbeResp = zeros(1, nTrl);
    for i=1:nTrl
        t = mod(i-1,actTrl)+1;

        trl(i).trlIdx = mat.results.trialOrder(1,t,b);

        mat.results.stimCenters(2,:,trl(i).trlIdx,b) = -mat.results.stimCenters(2,:,trl(i).trlIdx,b) + mat.dpP.winRect(4);

        % Extrai a atividade de movimento ocular
        trlEyeData = eyeData;
        trlEyeData.coord = eyeData.coord(:, find(eyeData.time ==trlLimAbs(1,i), 1):find(eyeData.time ==trlLimAbs(2,i), 1)) - screenCenter;
        if isempty(trlEyeData.coord) || sum(isnan(trlEyeData.coord(:)))/numel(trlEyeData.coord(:)) > 1/3
            warning('Trial sem dados suficientes! (cf. foragingTrlProps)')
            trlKeep = 0;
        else
            eyeMovs = microsacc_emd(trlEyeData, "fieldData", 'coord', "fieldPupil", 'pupil', "fieldFreq", 'fs', "dataUnits", 'pixel', ...
                "findEvents", {'sac', 'fix', 'blk'}, "screenRes", screenRes, "screenDist", screenDist, "screenWidth", screenWidth, ...
                "saveRawData", true, "saveFiltData", true);
    
            % Um trial não é mantido se anteriormente não era bom (o que
            % acontece quando aborta o índice por excesso de repetições),
            % ou porque não foi registrado movimento ocular necessário
            trlKeep = mat.results.trialOrder(2,t,b) && ~isempty(eyeMovs) && ~isempty(eyeMovs.fixations.lims) && ~isempty(eyeMovs.saccades.lims);
        end

        trl(i).probeSeen = mat.results.isSaccSeen(b,trl(i).trlIdx);
        if expIdx ~= 2
            trl(i).trlKeep = trlKeep && ~trl(i).probeSeen;
        else
            trl(i).trlKeep = trlKeep;
        end

        % Só entro aqui se a probe não tiver sido vista
        if trl(i).trlKeep

            % Todos os índices do trial
%             currTrlAllIdx = trlLimIdx(1,i):trlLimIdx(2,i);

            %% Gerais
            trl(i).subjNum = subj;
            trl(i).trlDur    = trlDur(i);
            trl(i).trlLimEvt = trlLimIdx(:,i); trl(i).trlLimClk = trlLimAbs(:,i);
%             trl(i).trlLimTime= getPhaseTime(eyeData.time, trl(i).trlLimClk, trl(i).trlLimClk);
            trl(i).blkIdx = b;
            trl(i).tgtOri = mat.results.targetOri(b); trl(i).nTgt   = mat.results.nTs(b,trl(i).trlIdx);
            trl(i).tgtIdx = find(mat.results.orientation(:,trl(i).trlIdx, b) == trl(i).tgtOri);
            trl(i).nonTgtIdx = setdiff(1:mat.tkP.nStims, trl(i).tgtIdx);
            trl(i).tgtPosPix = mat.results.stimCenters(:,trl(i).tgtIdx,trl(i).trlIdx,b) - screenCenter;
            trl(i).tgtPosDva = pixel_to_dva(trl(i).tgtPosPix, 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            trl(i).stmOri    = mat.results.orientation(:,trl(i).trlIdx, b); trl(i).nStm = mat.tkP.nStims;
            trl(i).stmPosPix = mat.results.stimCenters(:,:,trl(i).trlIdx,b) - screenCenter;
            trl(i).stmPosDva = pixel_to_dva(trl(i).stmPosPix, 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));

            modTimes = mat.results.modTimes(b, trl(i).trlIdx);
            % A partir dele, sei que RR aparece após modTimes fixações
            % distintas, e modTimes + 1 é a fixação em pré-sacádico. Por
            % outro lado, os forHistIdx retornam as quantidades
            % considerando o total de fixações em ROIs feitas,
            %% independentemente de serem repetidas ou não
            trl(i).nROIsBeforeUpdate = modTimes;

            %% Obtém os limites de todos os intervalos de fixação em estímulos:
            %% SEMPRE DELIMITADOS POR ESTÍMULOS OU FIXAÇÕES
            % stm: fixações em estímulo (i.e., ROI + presença de estímulo);
            % ROI: fixações nas regiões dos estímulos (não precisa haver estímulo)
            % fix: fixações individuais que constituem os ROIs
            [trlKeep, perPhase, perROI, perStm, limsTime, hasRepetition] = getFixStmPhaseLims1(messages, eventLimClk, stmMsgs, eyeMovs, trl(i).trlLimEvt, trl(i).trlLimClk, mat.prm, modTimes);

            if trlKeep
                trl(i).fixPerPhase = perPhase.fix;
                trl(i).ROIPerPhase = perPhase.ROI;
                trl(i).stmPerPhase = perPhase.stm;
    
                trl(i).fixPerROI = perROI.fix;
                trl(i).stmPerROI = perROI.stm;
    
                trl(i).fixPerStm   = perStm.fix;
    
                trl(i).phaseLimsTime = limsTime.phase/mat.prm.fs;
                trl(i).ROILimsTime   = limsTime.ROI  /mat.prm.fs;
                trl(i).stmLimsTime   = limsTime.stm  /mat.prm.fs;
                trl(i).fixLimsTime   = limsTime.fix  /mat.prm.fs;
            end

%             trl(i).phaseLimsIdx = phaseLimsTime; trl(i).stmLimsIdx = stmLimsTime; trl(i).fixLimsIdx = fixLimsTime;
            
            trl(i).P3EarlyStop = mat.results.isP3earlyStop(b,trl(i).trlIdx);
        end

        feedback = mat.results.trialFeedback{b,t};

        hasProbeResp(i) = ismember(0, feedback(2,:));
        trl(i).trlKeep = trlKeep && ~trl(i).probeSeen && hasProbeResp(i);
        if trl(i).trlKeep
             % Será a primeira sacada após o fim da última fixação que 
             % ocorre na fase 3, podendo ocorrer antes ou depois do 
             % início do alvo
            P3SaccIdx    = find(limsTime.sacc(1,:) >= limsTime.fix(2, find(perPhase.fix(3,:), 1, 'last')), 1, 'first');

            repQueue = stimsQueue{i}(1,:);
            phaseRepQueue = stimsQueue{i}(3,:);
            noRepQueue = unique(repQueue, 'stable');

            %% Detalhes de probe e sacada antecendente
            % Estímulo fixado durante a fase 3
            P3ForHistIdx = find(phaseRepQueue == 3, 1,'first'); P3StmIdx = repQueue(P3ForHistIdx);
            P3SaccLims   = eyeMovs.saccades.lims(:, P3SaccIdx);
            P3SaccAmp    = eyeMovs.saccades.amplitude(P3SaccIdx);
            P3SaccVPeak  = eyeMovs.saccades.v_peaks(P3SaccIdx);

            if limsTime.phase(3,2) > P3SaccLims(1)
                if limsTime.phase(3,2) > P3SaccLims(2)
                    trl(i).trlKeep = 0;
                    warning('Sacada terminando antes do ruído rosa!')
                else
                    trl(i).saccInterval = 0;
                end
            else
                trl(i).saccInterval = limsTime.phase(3,2) - P3SaccLims(1);
            end
            trl(i).saccLatency = limsTime.phase(3,1) - P3SaccLims(1);
            trl(i).saccAmpDva  = P3SaccAmp;
            trl(i).saccVelDvas = P3SaccVPeak;

            % Estímulo fixado pós modificação
            PMForHistIdx = find(phaseRepQueue == 4, 1,'first');
            PMStmIdx = repQueue(PMForHistIdx);
            PMFixIdx = find(limsTime.fix(1,:) >= P3SaccLims(2), 1, 'first');
            PMFixPosPix = dva_to_pixel(eyeMovs.fixations.pos(:, PMFixIdx), 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            PMTgtPosPix = mat.results.stimCenters(:, PMStmIdx, trl(i).trlIdx, b) - screenCenter;
            trl(i).saccAccPix  = PMFixPosPix - PMTgtPosPix;
            trl(i).probePosPix = PMTgtPosPix;
            trl(i).probePosFixPix = PMFixPosPix;
            trl(i).probeIdx = PMStmIdx;
            trl(i).probeFixDur = diff(double(limsTime.ROI(:, P3ForHistIdx+1)))/mat.prm.fs;

            trl(i).probeCat = any(trl(i).tgtIdx == PMStmIdx);
            trl(i).probeOri = trl(i).stmOri(PMStmIdx);
            trl(i).probeForHistIdx = PMForHistIdx;

            %% Detalhes de pré-probe e RR
            P3FixPosPix = dva_to_pixel(mean(eyeMovs.data.filt_lin(:, limsTime.ROI(1,P3ForHistIdx):limsTime.ROI(2,P3ForHistIdx)),2), 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            P3TgtPosPix = mat.results.stimCenters(:, P3StmIdx, trl(i).trlIdx, b) - screenCenter;
            trl(i).preProbeCat = any(trl(i).tgtIdx == P3StmIdx);
            trl(i).preProbeOri = trl(i).stmOri(P3StmIdx);

            % Atenção: a fixação no pré-probe tem duas fases: uma com
            % estímulos em tela, outra sem estímulos, isso se 
            trl(i).preProbeFixDur = diff(double(limsTime.ROI(:, P3ForHistIdx)))/mat.prm.fs;
            P3FixLimsTime = limsTime.ROI(:, perPhase.ROI(3,:));
            trl(i).P3FixLimsTime = P3FixLimsTime/mat.prm.fs;
            trl(i).P3ROIDurPerPhase   = [0; limsTime.phase(2,2)-P3FixLimsTime(1,1); min(limsTime.phase(3,2), P3FixLimsTime(2,end)) - limsTime.phase(3,1); max(0, P3FixLimsTime(2,end)-limsTime.phase(4,1))]/mat.prm.fs;
            trl(i).pinkNoiseDur = trl(i).P3ROIDurPerPhase(3);
            trl(i).preProbePosPix = P3TgtPosPix;
            trl(i).preProbePosFixPix = P3FixPosPix;
            trl(i).preProbeProbeDistDva = pixel_to_dva(vecnorm(trl(i).preProbePosPix - trl(i).probePosPix), 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            trl(i).preProbeIdx = P3StmIdx;
            trl(i).preProbeForHistIdx = P3ForHistIdx;

            nSStmIdx = feedback(1, feedback(2,:) == 1);
            nSStmPosPix = mat.results.stimCenters(:, nSStmIdx, trl(i).trlIdx, b) - screenCenter;
            trl(i).nSaccProbeOri = trl(i).stmOri(nSStmIdx);
            trl(i).nSaccProbeCat = any(trl(i).tgtIdx == nSStmIdx);
            trl(i).nSaccProbePosPix = nSStmPosPix;
            trl(i).preProbeNSaccProbeDistDva = pixel_to_dva(vecnorm(trl(i).preProbePosPix - trl(i).nSaccProbePosPix), 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            trl(i).nSaccProbeProbeDistDva    = pixel_to_dva(vecnorm(trl(i).nSaccProbePosPix - trl(i).probePosPix), 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            nSVec = trl(i).nSaccProbePosPix - trl(i).probePosPix;
            pPVec = trl(i).preProbePosPix - trl(i).probePosPix;
            trl(i).nSaccProbeProbeDistDeg = acosd(dot(nSVec, pPVec) / (norm(nSVec) * norm(pPVec)));
            trl(i).nSaccProbeIdx = nSStmIdx;

            forStmIdx = feedback(1, feedback(2,:) == -1);
            forStmPosPix = mat.results.stimCenters(:, forStmIdx, trl(i).trlIdx, b) - screenCenter;
            trl(i).forProbeOri = trl(i).stmOri(forStmIdx);
            trl(i).forProbeCat = any(trl(i).tgtIdx == forStmIdx);
            trl(i).forProbePosPix = forStmPosPix;

            %% Forrageamento inclui repetições, mas não o pré-sacádico
            % Por isso, sem repetições, terá comprimento modTimes
            trl(i).forHistIdx = repQueue(phaseRepQueue ~= 4);
            trl(i).forHistLen = numel(trl(i).forHistIdx);
            trl(i).forHistHasRep = hasRepetition;
            trl(i).forHistIdxNoRep = unique(trl(i).forHistIdx, 'stable');
            %% Duração das fixações em estímulo; 
            % por isso o primeiro não é o ROI MAS o último (durante P3) não
            % se limita à presença do estímulo (já que tem componentes em
            % P2, P3 e P4)
            trl(i).forHistFixDur = diff([limsTime.stm(:,1) limsTime.ROI(:, 2:end-1)])/mat.prm.fs;
            trl(i).forHistCat = ismember(trl(i).forHistIdx, trl(i).tgtIdx);
            trl(i).forHistOri = trl(i).stmOri(trl(i).forHistIdx);

            trl(i).forProbeFixDur  = trl(i).forHistFixDur(find(trl(i).forHistIdx == forStmIdx,1));
            %% O valor 0 indica o pré-probe. Como o forProbe nunca é o pré-probe, 
            % forProbeRecency >= 1
            trl(i).forProbeRecency = find(trl(i).forHistIdx == forStmIdx, 1, "last") - trl(i).forHistLen; % trl(i).forHistLen;
            trl(i).forProbeIdx = forStmIdx;
            trl(i).forProbeForHistIdx = find(trl(i).forHistIdx == forStmIdx, 1, 'last');


            ansObs = xor(ismember(feedback(1,:), trl(i).tgtIdx), ~feedback(3,:));
            trl(i).probeResp      = ansObs(feedback(2,:) == 0);
            trl(i).nSaccProbeResp = ansObs(feedback(2,:) == 1);
            trl(i).forProbeResp   = ansObs(feedback(2,:) == -1);
            trl(i).allResp = ansObs;

            trl(i).probeHit = feedback(3, feedback(2,:) == 0);
            trl(i).nSaccProbeHit = feedback(3, feedback(2,:) == 1);
            trl(i).forProbeHit   = feedback(3, feedback(2,:) == -1);
            trl(i).allHit        = feedback(3,:);
            
            trl(i).allProbesOrder = feedback(2,:);
            
        end
        if expIdx == 2 && t == mat.tkP.t1.tkP.nTrainTrials
            b = b+1;
        elseif expIdx ~= 2 && t == mat.tkP.nTrials
            b = b+1; 
        end
    end
    
end
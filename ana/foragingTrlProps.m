function [trl, eyeData, eventLimClk] = foragingTrlProps(mat, edf, sesStr, subj)

    expIdx = 4;
    N      = mat.tkP.nBlocks*mat.tkP.nTrials;
    trl    = getTrlStruct(N);

    %% Importa mensagens e eventos do arquivo edf
    messages = {edf.FEVENT(:).message}; messages(cellfun(@isempty, messages)) = {''};
%     events = {edf.FEVENT(:).codestring};
    eventLimClk  = [edf.FEVENT(:).sttime; edf.FEVENT(:).entime];

    % Extrai o que importa quanto à trajetória do olho
    Eye = mat.tkP.Eye;
    eyeData.coord = [edf.FSAMPLE.gx(Eye, :);edf.FSAMPLE.gy(Eye, :)];
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
    [trlLimIdx, nTrl, ~] = delBadTrl(mat, messages, sesLimIdx, trlLimIdx, nTrl);

    if nTrl ~=mat.tkP.nTrials*mat.tkP.nBlocks
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

    newCond = zeros(1, nTrl);
    for i=1:nTrl
        t = mod(i-1, mat.tkP.nTrials)+1;

        % Extrai a
        trlEyeData = eyeData;
        trlEyeData.coord = eyeData.coord(:, find(eyeData.time ==trlLimAbs(1,i), 1):find(eyeData.time ==trlLimAbs(2,i), 1)) - screenCenter;
        eyeMovs = microsacc_emd(trlEyeData, "fieldData", 'coord', "fieldPupil", 'pupil', "fieldFreq", 'fs', "dataUnits", 'pixel', ...
            "findEvents", {'sac', 'fix', 'blk'}, "screenRes", screenRes, "screenDist", screenDist, "screenWidth", screenWidth, ...
            "saveRawData", true, "saveFiltData", true);

        % Um trial não é mantido apenas se anteriormente não era bom (o
        % que acontece quando aborta o índice por excesso de repetições),
        % ou porque alguma das fases (p1 a pM) não apresenta a quantidade
        % necessária de estímulos
        trlKeep = mat.results.trialOrder(2,t,b);
        trl(i).trlKeep = trlKeep;
        if trl(i).trlKeep

            % Todos os índices do trial
%             currTrlAllIdx = trlLimIdx(1,i):trlLimIdx(2,i);

            %% Gerais
            trl(i).subjNum = subj;
            trl(i).trlDur = trlDur(i);
                                        trl(i).trlLimEvt = trlLimIdx(:,i); trl(i).trlLimClk = trlLimAbs(:,i);
            trl(i).trlLimTime= getPhaseTime(eyeData.time, trl(i).trlLimClk, trl(i).trlLimClk);
            trl(i).trlIdx = mat.results.trialOrder(1,t,b);
            trl(i).blkIdx = b;
            trl(i).tgtOri = mat.results.targetOri(b); trl(i).nTgt   = mat.results.nTs(b,trl(i).trlIdx);
            trl(i).tgtIdx = find(mat.results.orientation(:,trl(i).trlIdx, b) == trl(i).tgtOri);
            trl(i).tgtPosPix = mat.results.stimCenters(:,trl(i).tgtIdx,trl(i).trlIdx,b) - screenCenter;
            trl(i).tgtPosDva = pixel_to_dva(trl(i).tgtPosPix, 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            trl(i).stmOri = mat.results.orientation(:,trl(i).trlIdx, b); trl(i).nStm = mat.tkP.nStims;
            trl(i).stmPosPix = mat.results.stimCenters(:,:,trl(i).trlIdx,b);
            trl(i).stmPosDva = pixel_to_dva(trl(i).stmPosPix, 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));

            %% Obtém os limites de todos os intervalos de fixação em estímulos
            % bem como das fixações individuais que os constituem
            modTimes = mat.results.modTimes(b, trl(i).trlIdx);
            [trlKeep, fixLimsTime, stmLimsTime, stmPerPhase, phaseLimsTime, saccLimsTime, hasRep, stmLimsTimeRep] = getFixStmPhaseLims(messages, eventLimClk, stmMsgs, eyeMovs, trl(i).trlLimEvt, trl(i).trlLimClk, mat.prm, modTimes);
            trl(i).phaseLimsIdx = phaseLimsTime;
            trl(i).stmLimsIdx   = stmLimsTime;
            trl(i).fixLimsIdx   = fixLimsTime;
            trl(i).phaseLimsTime = phaseLimsTime/mat.prm.fs;
            trl(i).stmLimsTime = stmLimsTime    /mat.prm.fs;
            trl(i).fixLimsTime = fixLimsTime    /mat.prm.fs;
        end

        feedback = mat.results.trialFeedback{b,t};

        %% APÓS NOVA COLETA (modificar foragingGabors), REMOVER CONDIÇÃO EXTRA
        newCond(i) = ismember(0, feedback(2,:));
        trl(i).trlKeep = trlKeep && newCond(i);
        if trl(i).trlKeep
             % Será a primeira sacada após o início da fase 3, podendo
             % ocorrer antes ou depois do início do alvo
            P3SaccIdx    = find(saccLimsTime(1,:) >= phaseLimsTime(3,1), 1, 'first');


            P3StmIdx = stimsQueue{i}(1,(stmPerPhase == 3));
            P3SaccLims   = eyeMovs.saccades.lims(:, P3SaccIdx);
            P3SaccAmp    = eyeMovs.saccades.amplitude(P3SaccIdx);
            P3SaccVPeak  = eyeMovs.saccades.v_peaks(P3SaccIdx);

            trl(i).saccInt     = phaseLimsTime(3,2) - P3SaccLims(1);
            trl(i).saccInt1    = phaseLimsTime(3,1) - P3SaccLims(1);
            trl(i).saccAmpDva  = P3SaccAmp;
            trl(i).saccVelDvas = P3SaccVPeak;

            PMStmIdx = stimsQueue{i}(1,(stmPerPhase == 4));
            PMFixIdx = find(fixLimsTime(1,:) > P3SaccLims(2), 1, 'first');
            PMFixPosPix = dva_to_pixel(eyeMovs.fixations.pos(:, PMFixIdx), 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            PMTgtPosPix = mat.results.stimCenters(:, PMStmIdx, trl(i).trlIdx, b) - screenCenter;
            trl(i).saccAccPix  = PMFixPosPix - PMTgtPosPix;
            trl(i).probePosPix = PMTgtPosPix;
            trl(i).probePosFixPix = PMFixPosPix;

            trl(i).probeCat = any(trl(i).tgtIdx == PMStmIdx);
            trl(i).probeOri = trl(i).stmOri(PMStmIdx);
            trl(i).probeForHistIdx = modTimes;


            P3aux = find(stmPerPhase == 3);
            P3FixPosPix = dva_to_pixel(mean(eyeMovs.data.filt_lin(:, stmLimsTime(1,P3aux):stmLimsTime(2,P3aux)),2), 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));
            P3TgtPosPix = mat.results.stimCenters(:, P3StmIdx, trl(i).trlIdx, b) - screenCenter;
            trl(i).preProbeCat = any(trl(i).tgtIdx == P3StmIdx);
            trl(i).preProbeOri = trl(i).stmOri(P3StmIdx);
            trl(i).preProbeFixDur = diff(double(stmLimsTime(:, P3aux)))/mat.prm.fs;
            trl(i).pinkNoiseDur   = (stmLimsTime(2, P3aux) - phaseLimsTime(3,1))/mat.prm.fs;
            trl(i).preProbePosPix = P3TgtPosPix;
            trl(i).preProbePosFixPix = P3FixPosPix;
            trl(i).preProbeProbeDistDva = pixel_to_dva(vecnorm(trl(i).preProbePosPix - trl(i).probePosPix), 'dist', screenDist, 'width', screenWidth, 'res', screenRes(1));

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

            forStmIdx = feedback(1, feedback(2,:) == -1);
            forStmPosPix = mat.results.stimCenters(:, forStmIdx, trl(i).trlIdx, b) - screenCenter;
            trl(i).forProbeOri = trl(i).stmOri(forStmIdx);
            % Preciso calcular ainda
            trl(i).forProbeFixDur = [];
            trl(i).forProbeCat = any(trl(i).tgtIdx == forStmIdx);
            trl(i).forProbePosPix = forStmPosPix;

            trl(i).forHistLen = modTimes - 1;
            trl(i).forHistHasRep = hasRep;
            trl(i).forHistIdx = stimsQueue{i}(1,1:end-1);
            trl(i).forHistIdxNoRep = unique(trl(i).forHistIdx, 'stable');
            trl(i).forHistFixDur = diff(double(stmLimsTimeRep(:, 1:end-1)))/mat.prm.fs;
            trl(i).forHistCat = ismember(trl(i).forHistIdx, trl(i).tgtIdx);
            trl(i).forHistOri = trl(i).stmOri(trl(i).forHistIdx);

            ansObs = xor(ismember(mat.results.trialFeedback{b,t}(1,:), trl(i).tgtIdx), ~mat.results.trialFeedback{b,t}(3,:));
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
        if t == mat.tkP.nTrials, b = b+1; end
    end
    
end


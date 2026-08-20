function [trlKeep, perPhase, perROI, perStm, limsTime, hasRepetition] = getFixStmPhaseLims1(messages, eventLimClk, stmMsgs, eyeMovs, trlLimEvt, trlLimClk, prm, modTimes)
%     figure; plot(eyeMovs.data.filt_lin'); hold on; xline([phaseLimsTime(1,2); phaseLimsTime(2,2); phaseLimsTime(3,2); phaseLimsTime(4,2)], '-k', 'Linewidth', 1.5);
%     keepTrl = 1; saccLimsTime = [];
    % Dados os eventos de um trial e a posição ocular durante esse período,
    % identifica propriamente os inícios e fim de fixação em cada estímulo
    eventStartClk = eventLimClk(1, :); clear eventLimClk;
    
    % (1) Usa as mensagens do trial para delimitar grosseiramente quando
    %     ocorrem os estímulos. Assim, consegue distinguir fixações boas das ruins
          %% Lembre que manter fixação após RR conta como OLD
          % IN:  trlLimEvt;
          % OUT: trlMsgs e seus lims em Evt
    trlMsgs = messages(trlLimEvt(1):trlLimEvt(2));
    trlStmMsgsOn  = find(contains(trlMsgs, stmMsgs(1))) + trlLimEvt(1) - 1;
    trlStmMsgsOff = find(contains(trlMsgs, stmMsgs(2))) + trlLimEvt(1) - 1;
    trlStmMsgsLimsEvt = pairStartEnd(trlStmMsgsOn, trlStmMsgsOff, trlLimEvt);
    

    % (4) Encontra os limites das fases
    p1MsgsLimsEvt = [find(contains(trlMsgs, prm.msg.on.P1)) find(contains(trlMsgs, prm.msg.on.P2))] + trlLimEvt(1) - 1;  p1MsgsLimsTime= double(eventStartClk(p1MsgsLimsEvt) - trlLimClk(1)); p2MsgsLimsEvt = [find(contains(trlMsgs, prm.msg.on.P2)) find(contains(trlMsgs, prm.msg.on.P3))] + trlLimEvt(1) - 1;  p2MsgsLimsTime= double(eventStartClk(p2MsgsLimsEvt) - trlLimClk(1));
    p3MsgsLimsEvt = [find(contains(trlMsgs, prm.msg.on.P3)) find(contains(trlMsgs, prm.msg.off.P3))] + trlLimEvt(1) - 1; p3MsgsLimsTime= double(eventStartClk(p3MsgsLimsEvt) - trlLimClk(1)); pMMsgsLimsEvt = [find(contains(trlMsgs, prm.msg.on.PM)) find(contains(trlMsgs, prm.msg.off.PM))] + trlLimEvt(1) - 1; pMMsgsLimsTime= double(eventStartClk(pMMsgsLimsEvt) - trlLimClk(1));
    phaseLimsTime = [p1MsgsLimsTime; p2MsgsLimsTime; p3MsgsLimsTime; pMMsgsLimsTime];
    
    % (2) Obtém o tempo exato em que o olho entra e sai de cada estímulo,
    %     inclusive início e fim de sacadas, e obtém mensagens sobre a 
    %     utilidade de visita
          % IN:  trlMsgs e seus lims em Evt
          % OUT: lims e classificação dos stm

    trlStmMsgs = trlMsgs(trlStmMsgsLimsEvt - (trlLimEvt(1) - 1));
    stmLimsTime = eventStartClk(trlStmMsgsLimsEvt) - trlLimClk(1);
    badStmIdx = strcmp(trlStmMsgs(2,:), prm.msg.off.stm{2});               % Salva os ruins de fato
    shortStmIdx = diff(stmLimsTime) < 15;                                  % Os únicos curtos de fato a serem descartados  
                                                                           % poderão ser em PM, que não entram na fila mas aparecem no arquivo,
    auxShortStmIdx = shortStmIdx;                                                                       % e não como ruins...
    auxShortStmIdx(stmLimsTime(1,:) >= phaseLimsTime(3,2)) = 0;
    if any(auxShortStmIdx & ~badStmIdx)
        warning('Estímulos curtos detectados durante forrageamento!');
    end
    badStmIdx = badStmIdx | shortStmIdx;
    oldStmIdx = strcmp(trlStmMsgs(1,:), prm.msg.on.stm{2});
%     notGoodStmIdx = badStmIdx | oldStmIdx;   % Salva os ruins ou repetidos

    % (3) Descarta os stm ruins (pois não vou precisar deles nunca, não são
    %     contabilizados como fixação e em geral nem incluem uma fixação)
    stmLimsTime = stmLimsTime(:,~badStmIdx);
    oldStmLimsTime = stmLimsTime; %#ok<NASGU> 
    oldStmIdx   = oldStmIdx(~badStmIdx);
    
    

    if diff(phaseLimsTime(3,:))/prm.fs < prm.minP3Dur/2
        trlKeep = false;
        perPhase = []; perROI = []; perStm = []; limsTime = []; hasRepetition = [];
        warning('Fase 3 exibida curtíssima (cf. getFixStmPhaseLims1)')
        return
    end

    
    clear eventStartClk trlMsgs trlStmMsgsOn trlStmMsgsOff trlStmMsgsLimsEvt...
        p1MsgsLimsEvt p2MsgsLimsEvt p3MsgsLimsEvt pMMsgsLimsEvt...
        p1MsgsLimsTime p2MsgsLimsTime p3MsgsLimsTime pMMsgsLimsTime trlLimEvt...
        trlLimClk stmMsgs messages;

    % (4) Extrai relações entre fases, fixações, ROIs e stms
    fixLimsTime = eyeMovs.fixations.lims; 
    saccLimsTime = eyeMovs.saccades.lims;
    [trlKeep, perPhase, perROI, perStm, ROILimsTime, stmLimsTime, oldROIIdx] = updatePers(fixLimsTime, stmLimsTime, phaseLimsTime, oldStmIdx, prm);

    limsTime.fix   = fixLimsTime;
    limsTime.sacc  = saccLimsTime;
    limsTime.stm   = stmLimsTime;
    limsTime.ROI   = ROILimsTime;
    limsTime.phase = phaseLimsTime;
    hasRepetition = [];
    if ~trlKeep
        warning('Trial ruim: alguma fase está sem fixação! (cf. updatePers)')
        return; 
    end

    %% Formato do trial típico: modTimes em 2, exatamente 1 em 3 (o último de 2), 
    % ao menos 1 em 4
    hasRepetition = sum(perPhase.ROI(2,:)) > modTimes;
    badP2 = sum(perPhase.ROI(2,:)) < modTimes;
    badP3 = sum(perPhase.ROI(3,:)) ~= 1;
    badPM_1 = sum(~oldROIIdx & perPhase.ROI(4,:)) == 0;
    badPM_2 = sum(perPhase.fix(4,:)) == 0;
    
    if badP2 || badP3 || badPM_1 || badPM_2
        trlKeep = 0;
        if badP2
            warning('Trial ruim: menos estímulos em P2 que esperado (cf. getFixStmPhaseLims1)')
        elseif badP3
            warning('Trial ruim: P3 curtíssima (cf. getFixStmPhaseLims1)')
        elseif badPM_1
            warning('Trial ruim: sem estímulo bom pós-modificação (cf. getFixStmPhaseLims1)')
        elseif badPM_2
            warning('Trial ruim: sem fixação pós-modificação (cf. getFixStmPhaseLims1)')
        end
        return; 
    end
end


function [trlKeep, perPhase, perROI, perStm, ROILimsTime, stmLimsTime, oldROIIdx] = updatePers(fixLimsTime, stmLimsTime, phaseLimsTime, oldStmIdx, prm)
% Agora retorna 
%    perPhase: fix, ROI e stm;
%    perROI  : fix e stm;
%    perStm  : fix;
% ROILimsTime e stmLimsTime com limites temporais;
% oldROIIdx com apenas ROIs revisitados

    trlKeep = true;
    % Encontra todas as 
    %   (i) fixações que acontecem durante alguma fase;
    %  (ii) fixações que acontecem durantes algum stim;
    % (iii) stims que acontecem por fase;
    nFix      = size(fixLimsTime, 2);
    nPhases   = size(phaseLimsTime,1);
    stmLimsTime = double(stmLimsTime);

    %% --- Tratamento especial ao stim OLD da PM ---
    % Se for muito curto e não abranger nenhum fixação, deve estar sobre
    % uma sacada, logo quando o sujeito estava tirando os olhos do estímulo
    aux1StmLimsTime = stmLimsTime;
    aux1StmPerPhase = aux1StmLimsTime(1,:) < phaseLimsTime(:,2) & aux1StmLimsTime(2,:) > phaseLimsTime(:,1);
    aux1Old = find(aux1StmPerPhase(4,:) & oldStmIdx, 1, 'first');
    
    if ~isempty(aux1Old)
        cond1 = diff_events(biggest_intersection_events(aux1StmLimsTime(:, aux1Old), fixLimsTime)) <= 0;
        if (cond1 || diff(aux1StmLimsTime(:, aux1Old)) < 10)
            stmLimsTime(:, aux1Old) = [];
            oldStmIdx(aux1Old)   = [];
        end
    end
    clear aux1StmLimsTime aux1StmPerPhase aux1Old
    %% ---                  ---                  ---

    auxStmLimsTime = stmLimsTime;
    auxStmPerPhase = auxStmLimsTime(1,:) < phaseLimsTime(:,2) & auxStmLimsTime(2,:) > phaseLimsTime(:,1);

    %% Phase (I)
    perPhase.fix = false(nPhases,nFix);
    perPhase.fix(1,:) = fixLimsTime(1,:) > phaseLimsTime(1,2) - prm.minFixTime1*prm.fs & fixLimsTime(1,:) < phaseLimsTime(1,2);
    perPhase.fix(2:end,:) = fixLimsTime(2,:) > phaseLimsTime(2:end,1) & fixLimsTime(1,:) < phaseLimsTime(2:end,2);

    %% Stm

    % stms ajustados com base em fixações e fases: stms que duram mais de
    % uma fase são divididos
    auxStmIdx = find(sum(auxStmPerPhase) > 1);

    if isequal(find(auxStmPerPhase(:, auxStmIdx)),[2 3]')
        stmLimsTime = [stmLimsTime(:, 1:auxStmIdx-1) [stmLimsTime(1, auxStmIdx) phaseLimsTime(3,1); phaseLimsTime(2,2) phaseLimsTime(3,2)] stmLimsTime(:, auxStmIdx+1:end)];
    else
        if sum(auxStmPerPhase(3,:)) == 0
            warning('Trial ruim: fase 3 possivelmente curta, sem estímulo ou fixação ! (cf. updatePers)')
            trlKeep = false;
            perROI = []; perStm = []; ROILimsTime = []; oldROIIdx = [];
            return
        else
            warning('Trial duvidoso: estímulo ocorrendo em duas fases que não 2 e 3 (cf. updatePers)')
        end
    end

    nStms = size(stmLimsTime,2);
    perStm.fix = stmLimsTime(1,:)' < fixLimsTime(2,:) & stmLimsTime(2,:)' > fixLimsTime(1,:);

    for stm = 1:nStms
        currStmFix = perStm.fix(stm,:);

%        perPhase.fix(:, currStmFix)
        stmLimsTime(:,stm) = [max([stmLimsTime(1,stm) min(fixLimsTime(1,currStmFix))]);
                                 min([stmLimsTime(2,stm) max(fixLimsTime(2,currStmFix))])];
        if sum(currStmFix) == 0
            trlKeep = false; 
        end
    end

    if ~trlKeep
        perROI = []; perStm = []; ROILimsTime = []; oldROIIdx = [];
        return
    end

    %% ROI
        % O stmLimsTime atual mistura ROI e stm. Para virar puramente ROI, precisa
        % excluir o primeiro stm da P4 caso ele seja repetido (não será revisita 
        % pois estamos controlando para isso antes de entrar nessa função)
    % (1) ROI grosseiro
    
    auxIdx = find(auxStmPerPhase(end,:), 1); auxOldStm = oldStmIdx(auxIdx);
    
    if auxOldStm
        ROILimsTime = [auxStmLimsTime(:,1:auxIdx-2) [auxStmLimsTime(1,auxIdx-1); auxStmLimsTime(2,auxIdx)] auxStmLimsTime(:,auxIdx+1:end)];
        oldStmIdx(auxIdx) = [];
    else
        ROILimsTime = auxStmLimsTime;
    end
    oldROIIdx = oldStmIdx;
    clear auxStmPerPhase auxIdx auxOldStm oldStmIdx
    nROIs = size(ROILimsTime,2);

    % ROI ajustado com base nas fixações
    perROI.fix = ROILimsTime(1,:)' < fixLimsTime(2,:) & ROILimsTime(2,:)' > fixLimsTime(1,:);
    for roi = 1:nROIs
        currROIFix = perROI.fix(roi,:);
        ROILimsTime(:,roi) = [min([fixLimsTime(1,currROIFix)]);
                              max([fixLimsTime(2,currROIFix)])];
        if sum(currROIFix) == 0
            trlKeep = false; 
        end
    end
    perROI.stm = stmLimsTime(1,:) < ROILimsTime(2,:)' & stmLimsTime(2,:) > ROILimsTime(1,:)';

    %% Phase (II)
    perPhase.stm = stmLimsTime(1,:) <= phaseLimsTime(:,2) & stmLimsTime(2,:) >= phaseLimsTime(:,1);
    perPhase.ROI = ROILimsTime(1,:) <= phaseLimsTime(:,2) & ROILimsTime(2,:) >= phaseLimsTime(:,1);
end
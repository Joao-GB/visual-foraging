function [keepTrl, fixLimsTime, stmLimsTime, stmPerPhase, phaseLimsTime, fixPerPhase, saccLimsTime, hasRepP2, hasRepP4, stmLimsTimeRep] = getFixStmPhaseLims(messages, eventLimClk, stmMsgs, eyeMovs, trlLimEvt, trlLimClk, prm, modTimes)
%     figure; plot(eyeMovs.data.filt_lin'); hold on; xline([phaseLimsTime(1,2); phaseLimsTime(2,2); phaseLimsTime(3,2); phaseLimsTime(4,2)], '-k', 'Linewidth', 1.5);
    keepTrl = 1; saccLimsTime = [];
    % Dados os eventos de um trial e a posição ocular durante esse período,
    % identifica propriamente os inícios e fim de fixação em cada estímulo
    eventStartClk = eventLimClk(1, :); clear eventLimClk;
    
    %% Usa as mensagens do trial para delimitar grosseiramente quando
    % ocorrem os estímulos. Assim, consegue distinguir fixações boas das ruins
    trlMsgs = messages(trlLimEvt(1):trlLimEvt(2));
    trlStmMsgsOn  = find(contains(trlMsgs, stmMsgs(1))) + trlLimEvt(1) - 1;
    trlStmMsgsOff = find(contains(trlMsgs, stmMsgs(2))) + trlLimEvt(1) - 1;
    trlStmMsgsLimsEvt = pairStartEnd(trlStmMsgsOn, trlStmMsgsOff, trlLimEvt);
    
    % Obtém o tempo exato em que o olho entra e sai de cada ROI e as
    % mensagens sobre a utilidade de visita
    trlStmMsgs = trlMsgs(trlStmMsgsLimsEvt - (trlLimEvt(1) - 1));
    stmLimsTime = eventStartClk(trlStmMsgsLimsEvt) - trlLimClk(1);
    badStmOnlyIdx = strcmp(trlStmMsgs(2,:), prm.msg.off.stm{2});                                            % Salva os ruins de fato
    oldStmIdx = strcmp(trlStmMsgs(1,:), prm.msg.on.stm{2});
    badStmIdx = badStmOnlyIdx | oldStmIdx;   % Salva os ruins ou repetidos
    
    %% Encontra os limites das fases
    p1MsgsLimsEvt = [find(contains(trlMsgs, prm.msg.on.P1)) find(contains(trlMsgs, prm.msg.on.P2))] + trlLimEvt(1) - 1;  p1MsgsLimsTime= double(eventStartClk(p1MsgsLimsEvt) - trlLimClk(1)); p2MsgsLimsEvt = [find(contains(trlMsgs, prm.msg.on.P2)) find(contains(trlMsgs, prm.msg.on.P3))] + trlLimEvt(1) - 1;  p2MsgsLimsTime= double(eventStartClk(p2MsgsLimsEvt) - trlLimClk(1));
    p3MsgsLimsEvt = [find(contains(trlMsgs, prm.msg.on.P3)) find(contains(trlMsgs, prm.msg.off.P3))] + trlLimEvt(1) - 1; p3MsgsLimsTime= double(eventStartClk(p3MsgsLimsEvt) - trlLimClk(1)); pMMsgsLimsEvt = [find(contains(trlMsgs, prm.msg.on.PM)) find(contains(trlMsgs, prm.msg.off.PM))] + trlLimEvt(1) - 1; pMMsgsLimsTime= double(eventStartClk(pMMsgsLimsEvt) - trlLimClk(1));
    phaseLimsTime = [p1MsgsLimsTime; p2MsgsLimsTime; p3MsgsLimsTime; pMMsgsLimsTime];
    
    %% Com informação dos movimentos oculares, extrai fixações por fase,
    % bem como fixações por estímulo
    fixLimsTime = eyeMovs.fixations.lims;

    [fixPerPhase, fixPerStm, stmPerPhase] = updatePers(fixLimsTime, stmLimsTime, phaseLimsTime, prm);

    stmLimsTimeRep = stmLimsTime;

    %% Formato do trial típico: modTimes em 2, exatamente 1 em 3 (o último de 2), 
    % ao menos 1 em 4
    if sum(stmPerPhase(3,:)) > 1
        stmPerPhase(3,:) = stmPerPhase(3,:) & stmPerPhase(2,:);
    end
    hasRepP2 = sum(~badStmOnlyIdx & stmPerPhase(2,:)) ~= modTimes;
    hasRepP4 = any(stmPerPhase(4,:) & oldStmIdx);
    hasRepetition = hasRepP2 || hasRepP4;
    badP2 = sum(~badStmIdx & stmPerPhase(2,:)) ~= modTimes;
    badP3 = sum(~badStmIdx & stmPerPhase(3,:)) ~= 1;
    badPM_1 = sum(~badStmOnlyIdx & stmPerPhase(4,:)) == 0;
    badPM_2 = sum(fixPerPhase(4,:)) == 0;
    if badP2 || badP3 || badPM_1 || badPM_2
        keepTrl = 0;
        if badP2
            warning('Trial ruim: menos estímulos em P2 que esperado (cf. getFixStmPhaseLims)')
        elseif badP3
            warning('Trial ruim: P3 curtíssima (cf. getFixStmPhaseLims)')
        elseif badPM_1
            warning('Trial ruim: sem estímulo bom pós-modificação (cf. getFixStmPhaseLims)')
        elseif badPM_2
            warning('Trial ruim: sem fixação pós-modificação (cf. getFixStmPhaseLims)')
        end
        return; 
    end

    [keepTrl, stmLimsTimeAux, stmPerPhaseAux] = refineLims(eyeMovs, fixLimsTime, stmLimsTime, phaseLimsTime, fixPerPhase, fixPerStm, stmPerPhase, badStmIdx, prm);
    if hasRepetition
        [~, stmLimsTimeRep, ~] = refineLims(eyeMovs, fixLimsTime, stmLimsTime, phaseLimsTime, fixPerPhase, fixPerStm, stmPerPhase, badStmOnlyIdx, prm);
    else
        stmLimsTimeRep = stmLimsTimeAux;
    end
    stmLimsTime = stmLimsTimeAux;
    stmPerPhase = stmPerPhaseAux;

    if ~keepTrl
        warning('Trial ruim: P3 curta ou sem fixação nova em PM (cf. refineLims)')
        return; 
    end

    if sum(stmPerPhase<=3) ~= modTimes
        keepTrl = 0;
        warning('Trial ruim: menos estímulos em P2 que esperado (cf. getFixStmPhaseLims)')
        return; 
    end
    saccLimsTime = eyeMovs.saccades.lims;
end


function [fixPerPhase, fixPerStm, stmPerPhase] = updatePers(fixLimsTime, stmLimsTime, phaseLimsTime, prm)
    nSeenStms = size(stmLimsTime,2);
    nFix       = size(fixLimsTime, 2);
    stmPerPhase = zeros(4, nSeenStms);
    stmPerPhase(1,:) = stmLimsTime(1,:) < phaseLimsTime(1,2) & stmLimsTime(2,:) > phaseLimsTime(1,1);
    stmPerPhase(2,:) = stmLimsTime(1,:) < phaseLimsTime(2,2) & stmLimsTime(2,:) > phaseLimsTime(2,1);
    stmPerPhase(3,:) = stmLimsTime(1,:) < phaseLimsTime(3,2) & stmLimsTime(2,:) > phaseLimsTime(3,1);
    stmPerPhase(4,:) = stmLimsTime(1,:) < phaseLimsTime(4,2) & stmLimsTime(2,:) > phaseLimsTime(4,1);
    
    fixPerPhase = zeros(4, nFix);
    fixPerPhase(1,:) = fixLimsTime(1,:) > phaseLimsTime(1,2) - prm.minFixTime1*prm.fs & fixLimsTime(1,:) < phaseLimsTime(1,2);
    fixPerPhase(2,:) = fixLimsTime(2,:) > phaseLimsTime(2,1) & fixLimsTime(1,:) < phaseLimsTime(2,2);
    fixPerPhase(3,:) = fixLimsTime(2,:) > phaseLimsTime(3,1) & fixLimsTime(1,:) < phaseLimsTime(3,2); % Caso especial: a fixação deve terminar na fase 3, mas não precisa começar nela
    fixPerPhase(4,:) = fixLimsTime(2,:) > phaseLimsTime(4,1) & fixLimsTime(1,:) < phaseLimsTime(4,2);
    
    fixPerStm = zeros(nSeenStms, nFix);
    for i = 1:nSeenStms
        fixPerStm(i,:) = fixLimsTime(2,:) > stmLimsTime(1,i) & fixLimsTime(1,:) < stmLimsTime(2,i);
    end
end

function [keepTrl, stmLimsTime, stmPerPhase] = refineLims(eyeMovs, fixLimsTime, stmLimsTime, phaseLimsTime, fixPerPhase, fixPerStm, stmPerPhase, badStmIdx, prm)
    keepTrl = 1;
    % Se veio para cá, o trial é típico, então posso criar o stmPerPhase
    % desse jeito porco
    %% Refina limites para saber exatamente quando começou e terminou de ver 
    % cada estímulo, seja por movimento ocular ou por atualização da tela:
    % (1) a 1a fixação da P2 começa com o início da fase e termina no fim
    %     da última fixação que termina antes do fim do primeiro Stm;
    % (2) para os demais estímulos de completamente dentro de P2 fixações devem 
    %     começar depois e terminar antes do respectivo Stm
    % (3) Se houver transição que envolva fixação ruim em um estímulo mas em
    %     posição próxima à de uma fixação vizinha, considera que não houve
    %     movimento ocular, i.e., fixação faz parte do mesmo trlStm
    % (4) de P2 -> P3, se não há movimento ocular, o estímulo ficará fixado até
    %     desaparecer (fim de P3); se houver movimento ocular, quando a fixação
    %     for terminada, o estímulo desaparecerá e será o fim de P3. Em ambos
    %     os casos, o fim de P3 delimita o fim de um estímulo;
    % (5) de P3 -> PM, considero início de novo trlStm apenas se houver movi-
    %     mento ocular. Da mesma, será que fundo fixações espacialmente
    %     próximas?
    % (6) Os Stm inaceitáveis são complementares a reallyGoodStmIdx; os demais 
    %     talvez tenham sua utilidade

    %% P2
    trlStmLimsTime = [phaseLimsTime(2,1); fixLimsTime(2,find(fixPerStm(1, :), 1,'last'))];
    newStmPerPhase = [2]; %#ok<NBRAK2> 
    
    % Adiciona os demais da fase 2, i.e., nem primeiro nem último
    nSeenStimsP2 = sum(stmPerPhase(2,:));
    for i = 2:(nSeenStimsP2-1)
        if badStmIdx(i) == 1, continue; end
        trlStmLimsTime = [trlStmLimsTime [fixLimsTime(1,find(fixPerStm(i, :), 1,'first')); fixLimsTime(2,find(fixPerStm(i, :), 1,'last'))]]; %#ok<*AGROW> 
        newStmPerPhase = [newStmPerPhase 2];
    end
    %% P3 (último de P2)
    P3auxOn = fixLimsTime(1,find(fixPerStm(nSeenStimsP2, :), 1,'first')); P3auxOff = min(fixLimsTime(2,find(fixPerPhase(3, :), 1,'last')), phaseLimsTime(3,2));
    if isempty(P3auxOn) || isempty(P3auxOff)
        keepTrl = 0; return;
    end
    trlStmLimsTime = [trlStmLimsTime [P3auxOn; P3auxOff]];
    newStmPerPhase = [newStmPerPhase 3];

    %% P4: procuro pelo primeiro estímulo que não seja ruim, que ex
    SeenStimsP4Idx = find(stmPerPhase(4,:));

    count = 0;
    for i = SeenStimsP4Idx
        if badStmIdx(i) == 1, continue; end
        count = count+1;
        trlStmLimsTime = [trlStmLimsTime [max(fixLimsTime(1,find(fixPerStm(i, :), 1,'first')),phaseLimsTime(4,1)); min(fixLimsTime(2,find(fixPerStm(i, :), 1,'first')),phaseLimsTime(4,2))]];
        stmPerPhase = [newStmPerPhase 4];
    end
    if count == 0
        keepTrl = 0; return;
    end
    stmLimsTime = trlStmLimsTime;

%     aux = min(find(stmPerPhase(4, :), 1, 'first')+1, size(stmPerPhase, 2)); aux = find(fixPerStm(aux,:), 1, 'first'); if isempty(aux), aux = -Inf; else, aux = fixLimsTime(1, aux); end
%     stmLimsTime = [trlStmLimsTime [max(phaseLimsTime(4,1), aux); max(fixLimsTime(2,find(fixPerPhase(4, :), 1,'last')), phaseLimsTime(4,2))]];

%     nSeenStms = size(stmLimsTime,2);
%     nFix       = size(fixLimsTime, 2);
%     
% 
%     % Dadas duas fixações adjacentes, se dizem respeito à mesma posição
%     % espacial e ao menos um dos estímulos ao qual correspondem é ruim
% 
%     fixPerPos = zeros(1, nFix);
%     nPos = 1;
%     for i = 1:(nFix-1)
%         if any(sum(fixPerPhase(:, i:(i+1)),2)==2) && vecnorm(eyeMovs.fixations.pos(:,i)' - eyeMovs.fixations.pos(:,i+1)') <= prm.fixROIradius1_dva
%             fixPerPos(nPos, i:(i+1)) = 1;
%         else
%             fixPerPos(nPos, i) = 1;
%             fixPerPos = [fixPerPos; zeros(1, nFix)];
%             nPos = nPos + 1;
%             fixPerPos(nPos, i+1) = 1;
%         end
%     end
%     stmIdx = 1:nSeenStms;
%     stillBadIdx = badStmIdx;
%     wasBadIdx = badStmIdx;
%     for i = 1:(nSeenStms-1)
%         if sum(fixPerStm(:,i)) > 0 && sum(fixPerStm(:,i+1)) > 0
%             if all(fixPerPos(:,i) == fixPerPos(:,i+1)) && (wasBadIdx(i) || wasBadIdx(i+1))
%                 wasBadIdx(i+1) = 1;
%                 stillBadIdx(i:(i+1)) = 0;
%                 stmIdx(i+1)=stmIdx(i);
%             end
%         else
%             continue;
%         end
%     end
% 
%     s = find([true diff(stmIdx)~=0]);
%     e = [s(2:end)-1 numel(stmIdx)];
%     stmLimsTime = [stmLimsTime(1,s(~stillBadIdx(s))); stmLimsTime(2,e(~stillBadIdx(s)))];
end





    
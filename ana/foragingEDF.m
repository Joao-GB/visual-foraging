function [trlProps, eyeData, eventTimes] = foragingEDF(mat, edf, mode)
% Extrai todas as informações úteis do EDF por trial, 
    if nargin < 3, mode = 0; end

    % Mensagens do arquivo
    messages = {edf.FEVENT(:).message};
    messages(cellfun(@isempty, messages)) = {''};
    
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
        nTrl = size(trlLimIdx, 2);
    
        %% Encontra os trials descartados por pausa, timeout etc.
        badTrl = zeros(1, nTrl);
        for i=1:nTrl
            currTrlIdx = trlLimIdx(1,i):trlLimIdx(2,i);
            badCond = strcmp(messages(currTrlIdx), mat.prm.msg.err.trl{1}) | ... % Cobre efeitos de pausa
                strcmp(messages(currTrlIdx), mat.prm.msg.err.trl{2}) | ...
                strcmp(messages(currTrlIdx), mat.prm.msg.err.P1) | ...
                strcmp(messages(currTrlIdx), mat.prm.msg.err.P2);
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
    
        clear aux sesItMsg sesCpMsg auxIdx badblkIdx badBlk currBlkIdx blkOnIdx...
                goodTrl badTrl badCond currTrlIdx trlOnIdx trlOffIdx nBlk i
        stimQueue = reshape(mat.results.seenStimsQueue', 1, []);

        %% Extrai todas as fixações segundo as mensagens que escrevi no Eyelink
        % O limIdx é da forma: 
        % [início da fix; 
        % fim da fix;
        % duração da fix;         (exata, com base na 
        % % índice do estímulo visto segundo as matrizes de results;
        % visita a novo estímulo; (1 = sim)
        % fases em que ocorre;    (2 = P2, 3 = P3, 4 = PM, 5 = P2+P3, 7 = P3+P4)
        % término a força]        (1 = sim)

        stmLimIdx = cell(1, nTrl); P2LimIdx = zeros(2, nTrl); P3LimIdx = zeros(2, nTrl); PMLimIdx = zeros(2, nTrl);
        trlPerBlk = zeros(1, nTrl); fixDur = cell(1, nTrl); fixLims = cell(1, nTrl);

        events = {edf.FEVENT(:).codestring}; eventTimes = [edf.FEVENT(:).sttime; edf.FEVENT(:).entime];
        sampleTime = edf.FSAMPLE.time;

        for i=1:nTrl
            % Começo pelo fim pois ele me diz se a fixação é boa ou não
            trlPerBlk(i) = find(trlLimIdx(1,i)> blkIdx(1,:), 1);
            currTrlIdx = trlLimIdx(1,i):trlLimIdx(2,i);
            stmOffIdx = find(strcmp(messages(currTrlIdx), mat.prm.msg.off.stm{1})|strcmp(messages(currTrlIdx), mat.prm.msg.off.stm{3})|strcmp(messages(currTrlIdx), mat.prm.msg.off.stm{4})) + trlLimIdx(1, i) - 1;
            
            P2LimIdx(1:2,i) = [ find(strcmp(messages(currTrlIdx), mat.prm.msg.on.P2),1); find(strcmp(messages(currTrlIdx), mat.prm.msg.on.P3),1)] + trlLimIdx(1, i) - 1;
            P3LimIdx(1:2,i) = [ find(strcmp(messages(currTrlIdx), mat.prm.msg.on.P3),1); find(strcmp(messages(currTrlIdx), mat.prm.msg.off.P3),1)] + trlLimIdx(1, i) - 1;
            P3LimIdx(3, i)  = double(edf.FEVENT(P3LimIdx(2,i)).sttime - edf.FEVENT(P3LimIdx(1,i)).sttime)/1000;
    
            PMLimIdx(:,i) = [ find(strcmp(messages(currTrlIdx), mat.prm.msg.on.PM),1); find(strcmp(messages(currTrlIdx), mat.prm.msg.off.PM),1)] + trlLimIdx(1, i) - 1;
            
            auxIdx = [trlLimIdx(1,i) stmOffIdx];
            stmOnIdx = arrayfun(@(j) find(strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.on.stm{1})|strcmp(messages(auxIdx(j):auxIdx(j+1)), mat.prm.msg.on.stm{2}), 1,'last') + auxIdx(j) - 1,  1:length(stmOffIdx));

            % Encontra a primeira fixação a começar depois de cada stmOnIdx
            [fixDur{i}, fixLims{i}] = foragingGetFix(stmOnIdx, stmOffIdx, events, eventTimes);

            stmProperIdx = stimQueue{i}(1,:);
            newVisit = strcmp(messages(stmOnIdx), mat.prm.msg.on.stm{1});
            phaseOccur = 2*(stmOnIdx > P2LimIdx(1,i) & stmOnIdx < P2LimIdx(2,i)) + ...
                3*((stmOffIdx > P3LimIdx(1,i) & stmOffIdx < P3LimIdx(2,i))|(stmOnIdx > P3LimIdx(1,i) & stmOnIdx < P3LimIdx(2,i))) + ...
                4*(stmOffIdx > PMLimIdx(1,i) & stmOffIdx < PMLimIdx(2,i));
            
            forcedEnd = strcmp(messages(stmOffIdx), mat.prm.msg.off.stm{3})|strcmp(messages(stmOffIdx), mat.prm.msg.off.stm{4});
    
            stmLimIdx{i} = [stmOnIdx; 
                            stmOffIdx; 
                            fixDur{i};
                            stmProperIdx;
                            newVisit;
                            phaseOccur;
                            forcedEnd];
        end

        trlProps.ses  = sesLimIdx;
        trlProps.blk  = trlPerBlk; 
        trlProps.nTrl = nTrl;
        trlProps.trl  = trlLimIdx;
        trlProps.stm  = stmLimIdx;
        trlProps.fix  = fixLims;
%         trlProps.fixDur  = fixDur;
        trlProps.P2   = P2LimIdx;
        trlProps.P3   = P3LimIdx;
        trlProps.PM   = PMLimIdx;
        clear trlLimIdx trlPerBlk sesLimIdx P2LimIdx P3LimIdx ...
            PMLimIdx auxIdx currTrlIdx stmOnIdx stmOffIdx newVisit phaseOccur ...
            forcedEnd stimQueue P3stmIdx i PMminDur PMstmIdx

        % Extrai o que importa quanto à trajetória do olho
        Eye = mat.tkP.Eye;
        eyeData.coord = [edf.FSAMPLE.gx(Eye, :);edf.FSAMPLE.gy(Eye, :)];
        eyeData.time = edf.FSAMPLE.time;
    
        if mode
            %% Plota a trajetória ao redor da 3a fase
            isEyeEvent = strcmp(events,'ENDFIX') | strcmp(events,'ENDSACC');
            
%             eventTimesStart = [edf.FEVENT.sttime];
            eventTimesEnd   = [edf.FEVENT.entime];
            
            
        
            plotsPerFig = 10;
            nRows = 2; nCols = 5;
            figIdx = 0;
        
        
            % Demarcando os limites dessa fase e da fixação mais próxima em um estímulo
            for i = 1:nTrl
                if mod(i-1, plotsPerFig) == 0
                    figIdx = figIdx + 1;
                    figure('Name', sprintf('P3 traces – Fig %d', figIdx));
                end
            
                subplot(nRows, nCols, mod(i-1, plotsPerFig)+1); hold on
            
                % Limites temporais da 3a fase
                P3start = edf.FEVENT(trlProps.P3(1,i)).sttime;
                P3end   = edf.FEVENT(trlProps.P3(2,i)).sttime;
        
                PMstart = edf.FEVENT(trlProps.PM(1,i)).sttime;
                PMend   = edf.FEVENT(trlProps.PM(2,i)).sttime;
            
                % Encontra o índice do 1o evento terminado durante a 3a fase
                plotStIdx = find(isEyeEvent & eventTimesEnd > P3start, 1);
            
                % Mesma lógica para depois do fim da 3a fase
                plotEndIdx = find(isEyeEvent & eventTimesEnd > PMend, 1);
            
                % Limites temporais de cada subplot
                plotStTime  = edf.FEVENT(plotStIdx).sttime;
                plotEnTime = edf.FEVENT(plotEndIdx).entime;
            
                % Obtém em índices em vez de tempo
                s0 = find(sampleTime >= plotStTime, 1, 'first');
                s1 = find(sampleTime <= plotEnTime, 1, 'last');
            
                % Trajetória ocular
                x = edf.FSAMPLE.gx(Eye, s0:s1);
                y = edf.FSAMPLE.gy(Eye, s0:s1);
                t = double(sampleTime(s0:s1)) - double(P3start);    % O tempo é 0 no início da P3
                plot(t, x, 'k'); plot(t, y, 'k');
                title(sprintf('Trial %d', i));
            
                % Índices da fase 3 para demarcar seu início e fim no plot
        %         sP3start = double(find(sampleTime >= P3start, 1, 'first')) - double(P3start);
        %         sP3end   = double(find(sampleTime >= P3end,   1, 'first')) - double(P3start);
        %         sPMstart = double(find(sampleTime >= PMstart, 1, 'first')) - double(P3start);
        %         sPMend   = double(find(sampleTime >= PMend,   1, 'first')) - double(P3start);
        
                %% Para ver a trajetória dos olhos junto com a posição dos estímulo, basta adicionar pausa condicionada no trial
        %         stmCenters = mat.results.stimCenters(:, :,mat.results.trialOrder(1,i), 1);
        %         rad = mat.txP.gabor.size_px/2;
        %         for k = 1:size(stmCenters,2)
        %         rectangle('Position', [stmCenters(1,k)-rad, stmCenters(2, k)-rad, ...
        %         2*rad, 2*rad], ...
        %         'Curvature', [1 1], ...
        %         'EdgeColor', 'b', ...
        %         'LineWidth', 2);
        %         end
        %         set(gca, 'YDir', 'reverse');
        
                yl = ylim;
        
                plot([P3start, P3start] - P3start, yl, 'b-', 'LineWidth', 2);
                plot([P3end, P3end] - P3start, yl, 'b-', 'LineWidth', 2);
        
                plot([PMstart, PMstart] - P3start, yl, 'y-', 'LineWidth', 2);
                plot([PMend, PMend] - P3start, yl, 'y-', 'LineWidth', 2);
        % 
                % Demarco os limites das fixações que ocrrem na fase 3 ou PM em
                % vermelho
                for j=1:size(stmLimIdx{i}, 2)
                    if stmLimIdx{i}(6,j) == 4 || stmLimIdx{i}(6,j) == 5
                        goodStmStart = edf.FEVENT(stmLimIdx{i}(1,j)).sttime;
                        goodStmEnd   = edf.FEVENT(stmLimIdx{i}(2,j)).sttime;
                
                %         slStmStart = find(sampleTime >= goodStmStart, 1, 'first');
                %         slStmEnd   = find(sampleTime >= goodStmEnd, 1, 'first');
                
                        plot(double([goodStmStart, goodStmStart]) - double(P3start), yl, 'r--', 'LineWidth', 1);
                        plot(double([goodStmEnd, goodStmEnd]) - double(P3start), yl, 'r--', 'LineWidth', 1);
                    end
                end
                hold off
            end
            clear i j goodStmEnd goodStmStart P3end P3start PMend PMstart plotEnTime plotStTime plotEndIdx plotStIdx plotStTime...
                s0 s1 sP3end sP3start sPMstart sPMend t x y yl
            clear events eventTimesStart eventTimesEnd Eye figIdx isEyeEvent nCols nRows plotsPerFig sampleTime
        end
        
        
        %     trlOnIdx(logical(goodTrial))
        
                % Encontra os blocos bons da sessão experimental
                % Dos blocos bons, extrai apenas os trials bons
                % dos trials bons, extrai apenas as fixações completas que ocorreram
                % entre o início da fase 1 e o início da fase 3
        
            % Extrai os índices de início e fim dos trial bons
        %     blkOn = [find(cellfun(@(x) strcmp(x, mat.prm.msg.on.blk), messages))];
        %     goodTrlLims = 
        %     P3 = mov_props([find(cellfun(@(x) strcmp(x, 'TRIAL ONSET PHASE 3'), messages)); find(cellfun(@(x) strcmp(x, 'TRIAL OFFSET PHASE 3'), messages))]);
        %     evBetween = mov_props(trim_overlapping_events(P3.lims, [P3.lims(:)'; P3.lims(:)']));
        %     {edf.FEVENT(evBetween.idx).codestring}
           
        
    else
        disp('ERRO: não há sessões experimentais no arquivo!')
    end
        
        % 
        % auxMsg = {edf.FEVENT(2142:end).message};
        % auxMsg(cellfun(@isempty, auxMsg)) = {''};
        % 
        % trlOnIdx = find(cellfun(@(x) contains(x, 'TRIAL ONSET PHASE 2'), auxMsg));
        % trialOnsetIdx = trlOnIdx(1:5:end);
        % offsetIdx = find(cellfun(@(x) contains(x, 'TRIAL OFFSET'), auxMsg));
        % p3OffsetIdx = offsetIdx(1:3:end);
        % fixSearchLims = [trialOnsetIdx; p3OffsetIdx];
        % 
        % 
        % fixLenEDF = [];
        % codestring = {edf.FEVENT(:).codestring};
        % for i=1:length(fixSearchLims)
        %     auxIdx = find(cellfun(@(x) strcmp(x, 'ENDFIX'), codestring(fixSearchLims(1,i):fixSearchLims(2,i) ) ));
        %     auxIdx = auxIdx + fixSearchLims(1,i) - 1;
        %     % Recupera a duração de cada fixação
        %     fixLenEDF = [fixLenEDF [edf.FEVENT(auxIdx).entime]-[edf.FEVENT(auxIdx).sttime]]; %#ok<*AGROW> 
        % end
        % fixLenEDF = double(fixLenEDF);
        % figure; hist(fixLenEDF);
        % 
        % fixLenMine = cell2mat(reshape(results.seenStimsQueue, 1,[]));
        % fixLenMine = fixLenMine(2,:)*1000;
        % 
        % 
        % figure; hist(fixLenMine)
        
end
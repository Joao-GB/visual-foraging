% Versão só com 1 resposta do treino da tarefa: do início da fase 4 à linha
% ~2184. Copiei sem me importar com ifs e tudo mais, basta substituir se
% quiser recuperar

            % x. Desenha placeholders para os estímulos aleatoriamente, obedecendo
            %    a identificação deles como visitados, atual e não visitados
                        if mode == 1
                            Screen('Close', bg); clear auxWin bg;
                        end

                        % % Se houver fixado em algum estímulo em PM...
                        skipP4 = true;
                        if ~isempty(currIdx)
                            skipP4 = false;
                            fprintf('Há currIdx... ')
                            fprintf('então será perguntado sobre\n');
                        else
%                             nPre = nStimsToReport(1, idx, b); nPost = nStimsToReport(3, idx, b);
                            fprintf('Não há currIdx... ')
                            % if mode <= 2
                            %     if rand < prm.seenNotSeenRatio/(1+prm.seenNotSeenRatio)
                            %         fprintf('... então pré-s. vira visto\n')
                            %         nPre = nStimsToReport(1, idx, b) + nStimsToReport(2, idx, b);
                            %         nPost = nStimsToReport(3, idx, b);
                            %     else
                            %         fprintf('... então pré-s. vira não visto\n')
                            %         nPre = nStimsToReport(1, idx, b);
                            %         nPost = nStimsToReport(3, idx, b) + nStimsToReport(2, idx, b);
                            %     end
                            %     nStimsToReport(1, idx, b) = nPre;
                            %     nStimsToReport(2, idx, b) = 0;
                            %     nStimsToReport(3, idx, b) = nPost;
                            % 
                            %     if numel(seenIdx) < nPre
                            %         dif = nPre - numel(seenIdx);
                            %         nStimsToReport(1, idx, b) = numel(seenIdx);
                            %         nStimsToReport(3, idx, b) = nPost + dif;
                            %     end
                            %     if numel(notSeenIdx) < nPost
                            %         dif = nPost - numel(notSeenIdx);
                            %         nStimsToReport(3, idx, b) = numel(notSeenIdx);
                            %         nStimsToReport(1, idx, b) = nPre + dif;
                            %     end
                            % end
                        end
%                         disp(nStimsToReport(:, idx, b))
                        if sum(nStimsToReport(:, idx, b)) ~= 3
                            disp('ERRO')
                        end

                        % Para o treino de forrageamento, não faz sentido
                        % punir com o vermelhinho, até porque ainda está 
                        % sendo feito o ajuste dos parâmetros temporais
                        if mode <= 2, skipP4 = false; end

                        if skipP4
                            warningFlip(dpP.window, stimCenters(:, :, idx, b), resizeRect(dstRects, .25), currStim, txP.gabor.size_px, drP.allPW, drP.darkRed);
                        else
                            isTargetAnswer = nan(1,nStims); allColors2 = drP.allColors;
                            Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
                            % Os demais pontos que não o fixado antes da atual
                            nbhd = NeighborsOrder(stimCenters(:, :, idx, b), currStim, currIdx);
    
                            % O do passado não importa de onde eu pergunto,
                            % desde que não seja a última fixação nem o pós
                            % sacádico
                            seenAux = [];
                            if mode >= 3
                                isSaccSeen(b, idx) = any(ismember(currIdx, seenIdx));
                                auxIdx = setdiff(seenIdx, [currStim currIdx]);
                                seenAux = datasample(auxIdx, min(length(auxIdx), nStimsToReport(1, idx, b)), 'Replace', false);
                            end

                            currAux = []; 
                            if nStimsToReport(2, idx, b) == 1, currAux = currIdx; end

                            % Se nStimsTOReport somar 2 para essas duas
                            % condições, sempre haverá pelo menos um pra
                            % cada com o mínimo de vistos antes da
                            % atualização sendo 2
                            if isempty(seenAux) || isempty(currAux)
                                theAux = datasample(seenIdx, min(length(seenIdx), nStimsToReport(1, idx, b)+nStimsToReport(2, idx, b)), 'Replace', false);
                                seenAux = theAux(1:nStimsToReport(1, idx, b));
                                if mode >= 3
                                    currAux = theAux((end-nStimsToReport(2, idx, b)+1):end);
                                else
                                    currAux = [];
                                end
                            end
                            % O não visto eu pergunto o mais próximo de
                            % currStim que não seja currIdx e nem já visto. Mas
                            % se a pessoa nao mover o olho, como puniçao,
                            % perguntamos sobre os mais distantes
                            if isempty(currIdx)
                                nbhd = flip(nbhd);
                            end
    
                            notSeenNbhd = setdiff(nbhd, currIdx, 'stable');
                            notSeenNbhd = setdiff(notSeenNbhd, seenIdx, 'stable');
                            notSeenAux = notSeenNbhd(1:min(numel(notSeenNbhd),nStimsToReport(3, idx, b)));
    
                            % notSeenAux = datasample(auxNotSeenIdx, min(length(auxNotSeenIdx), nStimsToReport(3, idx, b)), 'Replace', false);
                            % Tanto seenAux como notSeenAux devem ser linhas
                            % para concatenar
                            fprintf('Vistos: '); fprintf(num2str(seenAux));
                            if mode >= 3
                                fprintf('\nPré-s: '); fprintf(num2str(currAux));
                                fprintf('\nNão vistos: '); fprintf(num2str(notSeenAux));
                            else
                                fprintf('\nPré-s: AUSENTE');
                                fprintf('\nNão vistos: AUSENTE');
                            end
                            fprintf('\n');

                            if mode >=3
                                orderToReportStimsCell = {seenAux, currAux, notSeenAux};
                                orderToReportStims = [orderToReportStimsCell{orderToReportSets(1, idx, b)} orderToReportStimsCell{orderToReportSets(2, idx, b)} orderToReportStimsCell{orderToReportSets(3, idx, b)}];
                                if numel(orderToReportStims) < 2 || numel(orderToReportStims) > 3
                                    fprintf('Algo de errado: tem que reportar: %d', numel(orderToReportStims));
                                end
                                orderRemapped = [];
                                for auxIdx = 1:3
                                    orderRemapped = [orderRemapped orderToReportMap(orderToReportSets(auxIdx, idx, b))*ones(1, length(orderToReportStimsCell{orderToReportSets(auxIdx, idx, b)}))]; %#ok<AGROW> 
                                end
                            else
                                orderToReportStims = seenAux;
                                orderRemapped = ones([1 numel(seenAux)]);
                            end
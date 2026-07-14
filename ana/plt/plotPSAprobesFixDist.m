function plotPSAprobesFixDist(trl, mat)
    % 1. Extrair as distâncias
    distSacc  = [trl.preProbeProbeDistDva];
    distNSacc = [trl.preProbeNSaccProbeDistDva];
    
    % Flexibidade quanto à quantidade de nBins (e.g. 3 para tercis)
    nBins = 4; 
    probs = linspace(0, 1, nBins + 1);
    edges = quantile(distSacc, probs); 
    edges = .5*round(edges/.5);
    edges(1) = 0; edges(end) = Inf;
    
    gridSize = nBins + 1;
    
    % Inicializa as células para a grade 4x4
    plotDataAcc  = cell(gridSize, gridSize);
    plotDataSens = cell(gridSize, gridSize);
    countData    = cell(gridSize, gridSize);
    
    calcDPrime = @(hit, total) norminv(min(max(hit/total, 0.001), 0.999)) * sqrt(2);
    
    % 2. Loop pela grade de bins (Core 3x3)
    for r = 1:nBins % Linhas: Distância do Sacádico
        for c = 1:nBins % Colunas: Distância do Não-sacádico
            mask = (distSacc >= edges(r) & distSacc < edges(r+1)) & ...
                   (distNSacc >= edges(c) & distNSacc < edges(c+1));
            
            [plotDataAcc{r,c}, plotDataSens{r,c}, countData{r,c}] = processCell(trl(mask), calcDPrime);
        end
    end
    
    % 3. Marginais das Colunas (Linha inferior: ignora distância do Sacádico)
    for c = 1:nBins
        mask = (distNSacc >= edges(c) & distNSacc < edges(c+1));
        [plotDataAcc{gridSize,c}, plotDataSens{gridSize,c}, countData{gridSize,c}] = processCell(trl(mask), calcDPrime);
    end
    
    % 4. Marginais das Linhas (Coluna direita: ignora distância do Não-sacádico)
    for r = 1:nBins
        mask = (distSacc >= edges(r) & distSacc < edges(r+1));
        [plotDataAcc{r,gridSize}, plotDataSens{r,gridSize}, countData{r,gridSize}] = processCell(trl(mask), calcDPrime);
    end
    
    % 5. Média Global (Canto inferior direito: todos os trials)
    [plotDataAcc{gridSize,gridSize}, plotDataSens{gridSize,gridSize}, countData{gridSize,gridSize}] = processCell(trl, calcDPrime);
    
    titleBase = 'Efeito da distância do pre-probe fixation';
    
    % Renderiza Acurácia
    renderPSAprobesFixDist(plotDataAcc, countData, edges, 'Acertos (%)', ...
        'PSA fixation distance - Accuracy', [titleBase, ' (Acurácia)'], mat.drP, false);
        
    % Renderiza Sensibilidade
    renderPSAprobesFixDist(plotDataSens, countData, edges, 'Sensibilidade (d'')', ...
        'PSA fixation distance - Sensitivity', [titleBase, ' (Sensibilidade)'], mat.drP, true);
end

% --- Função auxiliar local para evitar repetição de código ---
function [acc, sens, count] = processCell(subTrl, calcDPrime)
    if isempty(subTrl)
        acc = [0, 0]; sens = [0, 0]; count = 0;
        return;
    end
    
    [~, counts] = getPSAeffect(subTrl);
    total = counts(2,2); % Contagem de trials (igual para ambos no mesmo subTrl)
    count = total;
    
    if total > 0
        acc = [(counts(1,2)/total)*100, (counts(1,3)/total)*100];
        sens = [calcDPrime(counts(1,2), total), calcDPrime(counts(1,3), total)];
    else
        acc = [0, 0]; sens = [0, 0];
    end
end
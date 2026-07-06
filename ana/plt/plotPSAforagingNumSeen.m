function plotPSAforagingNumSeen(trl, drP)
    % plotPSAforagingNumSeen plots Sacc vs N-Sacc performance grouped by numSeen
    
    % 1. Extrair o número de itens vistos por trial
    numSeen = [trl.probeForHistIdx];
    uniqueSeen = 2:6; % Valores inteiros esperados de 2 a 6
    
    % Inicializar matriz para guardar as porcentagens de acertos
    % Linhas: cada valor de numSeen (2 a 6)
    % Colunas: 1 = Sacádico, 2 = Não-sacádico
    barMatrix = zeros(length(uniqueSeen), 2);
    
    % 2. Calcular o efeito PSA para cada subgrupo de numSeen
    for i = 1:length(uniqueSeen)
        currentVal = uniqueSeen(i);
        
        % Filtra a estrutura original de acordo com o numSeen atual
        mask = (numSeen == currentVal);
        
        if any(mask)
            subTrl = trl(mask);
            [~, subCounts] = getPSAeffect(subTrl);
            
            % Calcula as porcentagens [Sacádico, Não-sacádico]
            pct_s = (subCounts(1,2) / subCounts(2,2)) * 100;
            pct_n = (subCounts(1,3) / subCounts(2,3)) * 100;
            
            barMatrix(i, :) = [pct_s, pct_n];
        else
            barMatrix(i, :) = [0, 0]; % Se não houver dados para aquele N
        end
    end
    
    %% 3. Renderizar o Gráfico de Barras Agrupadas
    figure('Name', 'PSA num seen effect', 'Color', 'w', 'Position', [100 200 800 500]);
    
    % Plotando a matriz: agrupa por linhas automaticamente
    b = bar(barMatrix, 'grouped', 'EdgeColor', [0 0 0]);
    
    % Como barMatrix tem 2 colunas, 'b' terá 2 handles de série (b(1) e b(2))
    b(1).FaceColor = drP.darkBlue;   % Primeira barra de cada par: Sacádico
    b(2).FaceColor = drP.paleBrown;  % Segunda barra de cada par: Não-sacádico
    
    % Configurações estéticas dos eixos
    set(gca, 'TickDir', 'out', 'Box', 'off');
    grid on;
    
    % Nomear os agrupamentos no eixo X de 2 a 6
    set(gca, 'XTick', 1:length(uniqueSeen));
    xticklabels(string(uniqueSeen));
    
    xlabel('Número de estímulos vistos');
    ylabel('Acertos (%)');
    ylim([0 100]);
    title('Efeito pré-sacádico por quantidade de estímulos vistos');
    
    % Criar legenda simples para identificar o que é S e N dentro do cluster
    lgd = legend({'Sacádico', 'Não-sacádico'}, 'Location', 'northeast', 'Box', 'off');
end
% Larsson, Nystrom e Stridh (2013): 
% Ve
% 0) Mudar para degree
% 1) Preprocessamento: outliers/ruído, piscadas
%     Exclusão dos pontos que saem 1.5° além da tela
%     Piscadas afetam o entorno da piscada: **movimentos que lembram
%      sacadas verticais**: Início e fim da piscada como primeiro mínimo 
%      local antes e depois do fim da sacada. Máxima duração: 700 ms
%     One-sample spikes removidas com limiar de amplitude mínima,
%      velocidades logo antes devem ser menores que durante esse spike
%      (para não confundir com PSO)
% 2) Sacadas
%  -> Os limiares (do módulo) da aceleração para cada direção são obtidos
%      como \lambda* \sigma_i, i \in {x,y}, para uma primeira separação dos
%      movimentos bruscos de fixações e pursuits.
%     Constrói indicadora que assume valor 1 se ao menos uma das
%      acelerações supera o seu respectivo limiar, resultando em **intervalos
%      de sacadas aproximadas**
%     Duas potenciais sacadas que ocorrem muito próximas (separadas por
%      menos que t_min) são aglutinadas
%     Duração mínima de T ms
%  -> Onset e offset de cada sacada potencial são detectados pela identificação da
%      amostra de maior velocidade (de índice k) e vai para frente e para trás,
%      verificando se ao menos um de três critérios são satisfeitos por uma
%      quantidade mínima de amostras consecutivas
%      (a) Calcula a direção da trajetória a cada instante alpha(i) = 
%          arctan((x(i+1) - x(i))/(y(i+1) - y(i))). A direção média de uma 
%          sacada é definida como gamma = 1/3*(alpha(k - 1) + alpha(k) 
%          + alpha(k + 1)). Uma sacada começa ou termina se
%          abs( alpha(n) - gamma ) > delta para ao menos K = tk*fs amostras
%          consecutivas, com fs taxa de amostragem, e escolhe-se como a
%          fronteira o ponto mais próximo do instante k de vel máxima --
%          **desvio da direção média local**, talvez dê para capturar
%          tirando a média móvel 
%      (b) **Taxa de variação da direção** eps(i) = alpha(i) - alpha(i-1),
%          se abs( eps(i) ) > beta para N = tn*fs amostra consecutivas, e
%          escolhe-se como a fronteira o ponto desses que esteja mais distante
%          do instante k de vel máxima
%      (c) **Diferenças entre mudanças de direção (i.e., aceleração da direção)**
%          A. Como ni é calculado:
%             (1) pego os intervalos não sacádicos
%             (2) Removo o trending (fazem com PLR em blocos de 100 ms), obtendo resíduos
%             (3) Calculo a direção das diferenças dos resíduos, alpha'
%             (4) Calculo a distância entre direções consecutivas, epsilon'
%             (5) Pego todos epsilon' acima do limiar beta, resultando numa distribuição
%             (6) pego o 90-ésimo percentil dessa distribuição como ni
%          B. Como ni é usado:
%             (1) para os dados nos approximate saccade intervals, pega os instantes que satisfazem |epsilon| > beta, do critério anterior (não confundir com epsilon')
%             (2) Calcula as diferenças entre epsilons consecutivos, vou chamar de d_e
%             (3) Verifica se M pontos consecutivos de d_e são menores que ni
%             (4) Se for o caso, deve ser sacada 
%      (d) Velocidade do onset e do offset precisa ser 20% menor que a
%          velocidade máxima
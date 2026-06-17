function trace_heteroscedasticity(data, n_bins)
% Assume que data já está filtrado (sem valores gigantes não biológicos)

% Ordena os dados e os segmenta em n_bins bins, de modo a obter a variância
% para cada bin
    r = sqrt(data(1,:).^2+data(2,:).^2);

    r(isnan(r)) = [];
    % Divido os dados ordenados em n_bins de mesmo tamanho
    sorted_r = sort(r);

    L = length(sorted_r);
    edges = round(linspace(0, L, n_bins+1));       % Índices das bordas dos intervalos de igual tamanho

    
    half  = (edges(2) - edges(1))/2;        % Os plots serão feitos no meio do bin
    mids  = edges(1:end-1) + half;
    medians = sorted_r(round(mids));

    points = zeros([1, n_bins]);                 % Pontos de variânciaa serem plotados
    counts = zeros([1, n_bins]);                 % quantidade de pontos dos dados usados para obter aquela variância
    % Calculo a variância dos dados em cada intervalo
    for i = 1:n_bins
        counts(i) = edges(i+1) - edges(i);
        points(i) = std(sorted_r(edges(i)+1:edges(i+1))); % Somo 1 a edges(i) para evitar que use um mesmo dado em diferentes std's
    end

    figure;
    yyaxis left
    plot(medians, points);

    % Absolutamente dispensável, pois parece que os intervalos têm tamanhos
    % iguais a menos de 1 ou 2
%     yyaxis right
%     plot(medians, counts);


    yyaxis left
    title("Eye-trace Heteroscedasticity")
    xlabel("Radial eye position (dva)")
    ylabel("Standard deviation")
    
%     yyaxis right
%     ylabel("Number of points")

%     figure;
%     plot(data(1, 70000:110000))
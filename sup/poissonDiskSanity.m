mW = 52.4; sR = 1920; sD = 57;
minDist_dva = 4.5;
tolerance_dva = .5;

ROIparams = [960 540 480 270];
nStims = 10;
skipEllipse = true;

N = 10000;

minDist = dva2pix(sD, mW, sR, minDist_dva);
tolerance = dva2pix(sD, mW, sR, tolerance_dva);

nGenerated = zeros(1,N);    % Número de pontos por mapa
nnDist = [];                % Distância do vizinho mais próximo
diffDistClosest = [];       % Diferença entre os 2 vizinhos mais próximos de cada ponto

for i = 1:N

    [fixCenter, stimCenters, ~] = getStimLocations2(ROIparams, nStims-1, minDist, false, skipEllipse);
    
    allCenters = [fixCenter stimCenters];
    currN = size(allCenters,2);
    
    % skip if too few points (optional, or keep them)
    if currN < 3
        continue;
    end

    nGenerated(i) = currN;
    
    % --- pairwise distances ---
    D = pdist2(allCenters', allCenters');
    D(D==0) = inf;
    
    % --- nearest neighbor ---
    d1 = min(D, [], 2);
    nnDist = [nnDist; d1];
    
    % --- second nearest neighbor ---
    D_sorted = sort(D, 2);
    d2 = D_sorted(:,2);
    
    % --- difference ---
    diffDistClosest = [diffDistClosest; (d2 - d1)];
end

%% (1) success rate
meanN = mean(nGenerated);
stdN  = std(nGenerated);
fprintf('nGenerated stats: ');
fprintf('  mean = %.2f, ', meanN);
fprintf('  std  = %.2f\n', stdN);

figure;
histogram(nGenerated, 'BinMethod', 'integers');
xlabel('Number of generated points');
ylabel('Count');
title('Distribution of number of stimuli per screen');

%% (2) nearest neighbor stats
meanNN = mean(nnDist);
stdNN  = std(nnDist);

fprintf('NN distance: mean = %.2f px (%.2f * minDist), std = %.2f\n', ...
    meanNN, meanNN/minDist, stdNN);

figure;
histogram(nnDist, 30);
title('Nearest Neighbor Distance');
xlabel('Distance (px)');
ylabel('Count');

%% (3) difference between 1st and 2nd NN
meanDiff = mean(diffDistClosest);
stdDiff  = std(diffDistClosest);

fprintf('NN diff: mean = %.2f px, std = %.2f\n', meanDiff, stdDiff);

figure;
histogram(diffDistClosest, 30);
title('Difference: 2nd NN - 1st NN');
xlabel('Distance difference (px)');
ylabel('Count');
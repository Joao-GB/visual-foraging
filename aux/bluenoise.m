%% poisson_eye_tracking_layouts.m
% Generates point layouts that look evenly distributed without obvious rows,
% rings, or grids. Good starting point for eye-tracking stimuli.
%
% Output:
%   - one figure with layouts for 8, 9, 10, 11, 12 points
%   - a PDF saved in the current folder
%
% Works in a unit square [0,1] x [0,1].
% You can scale coordinates later to your screen/stimulus area.

clear; clc; close all;

rng(7);  % reproducible

counts = 8:12;

% Main generation settings
nCandidatesPerStep = 300;   % more = more even distributions
nRestarts          = 80;    % more = better final layout
edgeMargin         = 0.08;  % keep points away from borders
savePdf            = true;
pdfName            = 'poisson_eye_tracking_layouts.pdf';

allLayouts = cell(numel(counts),1);

figure('Color','w','Position',[100 100 1200 700]);

for i = 1:numel(counts)
    n = counts(i);

    pts = generateBestLayout(n, nCandidatesPerStep, nRestarts, edgeMargin);
    allLayouts{i} = pts;

    subplot(2,3,i);
    scatter(pts(:,1), pts(:,2), 140, 'filled', ...
        'MarkerFaceColor', [0 0.35 1], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    axis equal;
    xlim([0 1]);
    ylim([0 1]);
    set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
    title(sprintf('%d points', n), 'FontWeight', 'bold');
end

% Empty last subplot
subplot(2,3,6);
axis off;
text(0, 0.8, 'Generation method:', 'FontWeight', 'bold', 'FontSize', 12);
text(0, 0.6, sprintf(['Best-candidate blue-noise sampling\n' ...
                      '+ scoring for even spacing\n' ...
                      '+ weak penalty for visible alignment']), ...
                      'FontSize', 11);
xlim([0 1]); ylim([0 1]);

sgtitle('Poisson-disk / Blue-noise-like layouts for eye-tracking', ...
    'FontSize', 14, 'FontWeight', 'bold');

if savePdf
    exportgraphics(gcf, pdfName, 'ContentType', 'vector');
    fprintf('Saved PDF: %s\n', fullfile(pwd, pdfName));
end

% Print coordinates
for i = 1:numel(counts)
    fprintf('\n--- %d points ---\n', counts(i));
    disp(allLayouts{i});
end

%% -------- Local functions --------

function ptsBest = generateBestLayout(n, nCandidatesPerStep, nRestarts, edgeMargin)
% Generate many candidate layouts and keep the best-scoring one.

    bestScore = -Inf;
    ptsBest = [];

    for r = 1:nRestarts
        pts = bestCandidateSample(n, nCandidatesPerStep, edgeMargin);
        score = layoutScore(pts);

        if score > bestScore
            bestScore = score;
            ptsBest = pts;
        end
    end
end

function pts = bestCandidateSample(n, nCandidatesPerStep, edgeMargin)
% Mitchell-style best candidate sampling:
% Each new point is chosen from many random candidates as the one farthest
% from existing points.

    pts = zeros(n,2);

    % First point
    pts(1,:) = edgeMargin + (1 - 2*edgeMargin) * rand(1,2);

    for k = 2:n
        candidates = edgeMargin + (1 - 2*edgeMargin) * rand(nCandidatesPerStep, 2);

        bestIdx = 1;
        bestVal = -Inf;

        for c = 1:nCandidatesPerStep
            p = candidates(c,:);

            % Distance to existing points
            d = sqrt(sum((pts(1:k-1,:) - p).^2, 2));
            minDist = min(d);

            % Small edge preference toward center to avoid border crowding
            edgeClearance = min([p(1), p(2), 1-p(1), 1-p(2)]);

            % Candidate value: mostly maximize spacing, slightly prefer not too close to edge
            val = minDist + 0.15 * edgeClearance;

            if val > bestVal
                bestVal = val;
                bestIdx = c;
            end
        end

        pts(k,:) = candidates(bestIdx,:);
    end
end

function s = layoutScore(pts)
% Higher score = more even-looking distribution with less obvious structure.
%
% Terms:
%   + large minimum distance
%   + low variation in nearest-neighbor distance
%   + mild penalty for collinear/aligned triples
%   + mild penalty for strong crowding near edges

    D = pairwiseDistances(pts);
    n = size(pts,1);

    % nearest-neighbor distances
    nn = zeros(n,1);
    for i = 1:n
        di = D(i,:);
        di(i) = Inf;
        nn(i) = min(di);
    end

    minNN   = min(nn);
    meanNN  = mean(nn);
    stdNN   = std(nn);
    cvNN    = stdNN / meanNN;

    % Edge penalty: discourage too many points near borders
    edgeDist = min([pts(:,1), pts(:,2), 1-pts(:,1), 1-pts(:,2)], [], 2);
    edgePenalty = mean(1 ./ (edgeDist + 0.02));

    % Alignment penalty: penalize nearly collinear triples
    alignPenalty = alignmentPenalty(pts);

    % Final score
    s = ...
        4.0 * minNN ...          % spread out
      - 1.8 * cvNN ...           % uniform nearest-neighbor distances
      - 0.015 * edgePenalty ...  % avoid too much border hugging
      - 0.8 * alignPenalty;      % avoid obvious rows
end

function D = pairwiseDistances(pts)
% Full pairwise Euclidean distance matrix without toolboxes.

    x = pts(:,1);
    y = pts(:,2);

    dx = x - x.';
    dy = y - y.';
    D = sqrt(dx.^2 + dy.^2);
end

function p = alignmentPenalty(pts)
% Penalize triples that are too close to being collinear.
% This helps reduce obvious row-like structure.

    n = size(pts,1);
    p = 0;

    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n
                a = pts(j,:) - pts(i,:);
                b = pts(k,:) - pts(i,:);

                na = norm(a);
                nb = norm(b);

                if na < 1e-12 || nb < 1e-12
                    continue;
                end

                % sine of angle between vectors
                sinTheta = abs(a(1)*b(2) - a(2)*b(1)) / (na*nb);

                % penalize if nearly collinear
                if sinTheta < 0.12
                    p = p + (0.12 - sinTheta);
                end
            end
        end
    end
end
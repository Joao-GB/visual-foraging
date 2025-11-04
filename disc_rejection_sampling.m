function points = disc_rejection_sampling(a, b, N, c, ROIshape, ROIparams, avoidPoints, plot_discs)
% Sorteia uniformemente, em [a1,b2]x[b1,b2], N discos de raio c, de modo que
% eles não se sobreponham. Usa um algoritmo de aceitação/rejeição simples.
    rng('shuffle');

    if N == 0, points = []; return; end

    if length(a) < 2, a = [0 a]; end
    if length(b) < 2, b = [0 b]; end

    if nargin < 5, ROIshape = 'square'; end
    if nargin < 6, ROIparams  = [];     end
    if nargin < 7, avoidPoints =[];     end
    if nargin < 8, plot_discs = false;  end

    if isempty(avoidPoints), thereArePointsToAvoid = false; 
    else, thereArePointsToAvoid = true; l = size(avoidPoints,2); end

    changeShape = false;
    if strcmp(ROIshape, 'circle') || strcmp(ROIshape, 'ellipse')
        changeShape = true;
        if length(ROIparams) < 1, ROIparams = [(a(2)+a(1))/2, (b(2)+b(1))/2, (a(2)-a(1))/2, (b(2)-b(1))/2];
        elseif length(ROIparams) < 4, ROIparams(4) = ROIparams(3);
        end
        c1 = ROIparams(1); c2 = ROIparams(2);
        s1 = ROIparams(3); s2 = ROIparams(4); 
        ellipseCurve = @(x,y) (x-c1).^2./s1^2+(y-c2).^2./s2^2 <= 1;
    end

    max_attempts = 1000 * N;
    points = zeros(N, 2);
    curr = 0;
    
    for i = 1:max_attempts
        if curr >= N
            break;
        end
        
        candidate = [(a(2)-a(1))*rand()+a(1), (b(2)-b(1))*rand()+b(1)];
        
        % Verifica a distância entre o ponto e os já sorteados
        valid = true;
        if changeShape
            valid = ellipseCurve(candidate(1),candidate(2));
        end
        if valid == true && curr > 0
            if min(vecnorm(points(1:curr,:)' - repmat(candidate, curr, 1)')) < c
                valid = false;
            end
        end
        if valid == true && thereArePointsToAvoid
            if min(vecnorm(avoidPoints - repmat(candidate, l, 1)')) < c
                valid = false;
            end
        end

        
        if valid, points(curr+1, :) = candidate; curr = curr + 1; end
    end
    
    if curr < N
        warning('Only generated %d points out of %d requested', curr, N);
    end

    if plot_discs
        plot_points_with_discs(points, a, b, c, changeShape, ROIparams);
    end

    points = points';

end

function plot_points_with_discs(points, a, b, c, Rs, Rp)
    % PLOT_POINTS_WITH_DISCS Plot points with discs of radius c/2
    
    figure;
    hold on;
    
    % Plot the domain
    if Rs == true
        curvature = [1 1]; 
        pos = [Rp(1)-Rp(3), Rp(2)-Rp(4), 2*Rp(3), 2*Rp(4)];
    else
        curvature = 0; 
    pos = [a(1), b(1), a(2), b(2)];
    end
    rectangle('Position', pos, 'Curvature', curvature, 'EdgeColor', 'k', 'LineWidth', 2);
    
    % Plot points
    scatter(points(:,1), points(:,2), 50, 'r', 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
    
    % Plot discs (radius = c/2)
    theta = linspace(0, 2*pi, 100);
    for i = 1:size(points, 1)
        x_disc = points(i,1) + (c/2) * cos(theta);
        y_disc = points(i,2) + (c/2) * sin(theta);
        plot(x_disc, y_disc, 'b--', 'LineWidth', 1, 'HandleVisibility','off');
    end
    
    % Formatting
    axis equal;
    xlim([a(1)-c, a(2)+c]);
    ylim([b(1)-c, b(2)+c]);
    grid on;
    title(sprintf('%d Points with Discs (radius = %.2f)', size(points,1), c));
    xlabel('X');
    ylabel('Y');
    
    % Add legend
    plot(nan, nan, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Points');
    plot(nan, nan, 'b--', 'DisplayName', sprintf('Discs (r=%.2f)', c));
%     legend('show');
    
    hold off;
end
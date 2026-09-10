function fixPointHeatmap(ePos, oPos, sgn, rsp, stmDiam, maxDiam, plotMode)
% ePos : [2 x N] centros dos estímulos (x,y)
% oPos: ou
%        (A) [2 x N] centros das fixações, ou
%        (B) struct com campos:
%            oPos.center : [2 x N]
%            oPos.spread : [2 x 2*N]  (cada par de colunas é covariância 2x2)
% sgn  : [1 X N] 1 = alvo, 0 = distrator
% rsp  : [1 x N] 1 = correto, 0 = incorreto
% stmDiam : diâmetro do estímulo
% maxDiam : máximo diâmetro de uma fixação

if nargin < 7 || isempty(plotMode)
    plotMode = 'points';
end

% Fixações recentralizadas com base na posição do estímulo
if isstruct(oPos)
    centers = oPos.center;
else
    centers = oPos;
end

relPos = centers - ePos;

x = relPos(1,:);
y = relPos(2,:);

figure; hold on; axis equal;
xlabel('Posição horizontal em relação ao centro');
ylabel('Posição vertical em relação ao centro');
title('Distribuição de fixações em relação ao centro dos estímulos');

% Círculos de referência
theta = linspace(0, 2*pi, 400);

plot((stmDiam/2)*cos(theta), (stmDiam/2)*sin(theta), 'b-', 'LineWidth', 1.2);
plot((maxDiam/2)*cos(theta), (maxDiam/2)*sin(theta), 'b--', 'LineWidth', 1.0);

sig = sgn==1; noi = sgn==0;
cor = rsp==1; err = rsp==0;
if strcmpi(plotMode,'heatmap')

    % ---- recentered positions already computed ----
    % x, y vectors exist

    % define subplot structure:
    % Row 1 = correct
    % Row 2 = incorrect
    % Col 1 = signal
    % Col 2 = noise

    figure;

    condMatrix = {
        sig & cor, 'Alvo correto';
        sig & err, 'Alvo incorreto';
        noi & cor, 'Distrator correto';
        noi & err, 'Distrator incorreto'
    };

    for k = 1:4

        subplot(2,2,k); hold on; axis equal;

        idx = condMatrix{k,1};

        if any(idx)

            % 2D histogram density
            nbins = 40;

            edges = linspace(-maxDiam/2, maxDiam/2, nbins);

            N = histcounts2(x(idx), y(idx), edges, edges);

            imagesc(edges, edges, N');
            set(gca,'YDir','normal');
            colormap hot;
            colorbar;

        end

        % ---- reference circles ----
        theta = linspace(0,2*pi,400);

        plot((stmDiam/2)*cos(theta), ...
             (stmDiam/2)*sin(theta), ...
             'b-', 'LineWidth',1.2);

        plot((maxDiam/2)*cos(theta), ...
             (maxDiam/2)*sin(theta), ...
             'b--', 'LineWidth',1.0);

        title(condMatrix{k,2});
        xlabel('Horizontal');
        ylabel('Vertical');

        xlim([-maxDiam/2 maxDiam/2]);
        ylim([-maxDiam/2 maxDiam/2]);

        grid on;
    end

    sgtitle('Mapas de densidade das fixações');

    return;  % stop here, do not execute point mode
end

% Para plotar as elipses
drawEllipse = @(mu, C, style) ...
    plotEllipse(mu, C, style);

if ~isstruct(oPos)
    h1 = plot(nan,nan,'o','Color','k','MarkerFaceColor','k','LineWidth', 1.2); %#ok<*NASGU> 
    h2 = plot(nan,nan, 'o','Color','k','MarkerFaceColor','w','LineWidth', 1.2);
    h3 = plot(nan,nan, 'o', 'Color','r','MarkerFaceColor','r','LineWidth', 1.2);
    h4 = plot(nan,nan, 'o', 'Color','r','MarkerFaceColor','w','LineWidth', 1.2);
    plot(x(sig & cor), y(sig & cor),'o','Color','k','MarkerFaceColor','k','LineWidth', 1.2);
    plot(x(sig & err), y(sig & err), 'o','Color','k','MarkerFaceColor','w','LineWidth', 1.2);
    plot(x(noi & cor), y(noi & cor), 'o', 'Color','r','MarkerFaceColor','r','LineWidth', 1.2);
    plot(x(noi & err), y(noi & err), 'o', 'Color','r','MarkerFaceColor','w','LineWidth', 1.2);

else
    h1 = plot(nan,nan,'Color','k','LineStyle','-','LineWidth', 1.2);%'Color','k','MarkerFaceColor','k'
    h2 = plot(nan,nan, 'Color','k','LineStyle','--','LineWidth', 1.2);%'Color','k','MarkerFaceColor','w'
    h3 = plot(nan,nan, 'Color','r','LineStyle','-','LineWidth', 1.2);%'Color','r','MarkerFaceColor','r'
    h4 = plot(nan,nan, 'Color','r','LineStyle','--','LineWidth', 1.2);%'Color','r','MarkerFaceColor','w'
    for i = 1:size(centers,2)

        mu = relPos(:,i);
        C  = oPos.spread(:, (2*i-1):(2*i)); % 2x2 covariance

        if sig(i) && cor(i)
            drawEllipse(mu, C, struct('Color','k','LineStyle','-'));
        elseif sig(i) && err(i)
            drawEllipse(mu, C, struct('Color','k','LineStyle','--'));
        elseif noi(i) && cor(i)
            drawEllipse(mu, C, struct('Color','r','LineStyle','-'));
        else
            drawEllipse(mu, C, struct('Color','r','LineStyle','--'));
        end
    end
end

legend({'Diâmetro estímulo', 'Tolerância', ...
        'Alvo correto', 'Alvo incorreto', ...
        'Distr. correto', 'Distr. incorreto'}, ...
        'Location','bestoutside');

grid on;
hold off;
end

function h = plotEllipse(mu, C, style)

% mu : 2x1 mean
% C  : 2x2 covariance
% style : struct with fields Color and LineStyle

[U,S,~] = svd(C);

t = linspace(0,2*pi,200);
circle = [cos(t); sin(t)];

ellipse = U * sqrt(S) * circle;   % scale unit circle by covariance
ellipse = ellipse + mu;           % shift to center

h = plot(ellipse(1,:), ellipse(2,:), ...
         'Color', style.Color, ...
         'LineStyle', style.LineStyle, ...
         'LineWidth', 1.2);
end
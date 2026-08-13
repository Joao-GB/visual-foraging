% Dado um sistema dinâmico definido em um espaço (no mínimo?) 2D, plota sua
% trajetória ao longo do tempo, dando para visualizar apenas a atividade
% recente ou total...
% Dá para melhorar muito:
% 1) Corrigir o aspecto/adaptabilidade das barras e textos ao tamanho do plot;
% 2) Adicionar uma forma geométrica padrão (elipses, retângulos, triângulos,
%    com dimensões a especificar) ou particular (dada em coordenadas) --
%    -- e.g. um **atrator**

% Disclaimer: as funções de botão e slider foram primeiro propostas pelo
%             ChatGPT, e modificações foram feitas em cima.
function interactive_trajectory(data, fs, varargin)
    % data: array 2D em que a primeira linha são as coordenadas em x e a
    % segunda em y
    % fs: taxa amostral com que os dados foram coletados
    
    x = data(1, :);
    y = data(2, :);
    
    options = set_options(varargin);
    x_w = max(abs([max(x),min(x)]));
    y_w = max(abs([max(y),min(y)]));

    % Altera o do plot principal conforme os dados
    xpy = x_w/y_w;
    T = length(x)/fs;
    h = 600;
    w = h*xpy;

    % Cria a figura inicial
    fig = figure('Position', [100 100 w h]);
    ax = axes(fig);
    plot(ax, x(1), y(1), 'o', 'MarkerFaceColor', options.color);  % Initial point
    
    for_every_plot(x(1));
    
    % Slider
    slider = uicontrol('Style', 'slider', ...
        'Position', [100 10 .6*w 20], ...
        'Min', 1, 'Max', length(x), ...
        'Value', 1, ...                                 % Valor inicial
        'SliderStep', [2/length(x) 100/length(x)], ...  % Step quando clico na seta e quando no espaço
        'Callback', @updatePlot);
    
    % Play/Pause Button
    playButton = uicontrol('Style', 'pushbutton', ...
        'Position', [.6*w+100 10 .1*w 20], ...
        'String', 'Play', ...
        'Callback', @togglePlay);
    
    % Add label for time
    timeLabel = uicontrol('Style', 'text','FontSize', 11, ...
        'Position', [.7*w+50 80 .15*w 20], ...
        'String', sprintf('Time: %.3f of %.1f s', 0, T));
    
    
    % Global variables for animation control
    isPlaying = false;
    currentIdx = 1;

    % Slider callback
    function updatePlot(src, ~)
        currentIdx = round(src.Value);
        t = currentIdx/fs;
        m = max(1, currentIdx-options.thr);

        if strcmp(options.past, "keep")
            plot(ax, x(1:currentIdx), y(1:currentIdx), 'b-', 'LineWidth', 1.5);
        else
            if strcmp(options.past, "fade")
                plot(ax, x(m:currentIdx), y(m:currentIdx), 'b-', 'LineWidth', 1.5);
            elseif strcmp(options.past, "lighter")
                plot(ax, x(1:m), y(1:m), 'b-', 'LineWidth', 1.5, 'Color', [0.357 0.812 0.957]);
                hold on;
%                 disp(m)
                plot(ax, x(m:currentIdx), y(m:currentIdx), 'b-', 'LineWidth', 1.5);
            end
            
            
        
%             plot(ax, x(:currentIdx), y(max(1, currentIdx-options.thr):currentIdx), 'b-', 'LineWidth', 1.5);
        end
        hold(ax, 'on');
        plot(ax, x(currentIdx), y(currentIdx), 'o', 'MarkerFaceColor', options.color);
        for_every_plot(x(currentIdx));
        hold(ax, 'off');
        timeLabel.String = sprintf('Time: %.3f of %.1f s', t, T);
    end
    
    % Play button callback
    function togglePlay(~, ~)
        isPlaying = ~isPlaying;  % Toggle play state
        if isPlaying
            playButton.String = 'Pause';
            for k = currentIdx:options.play_step:length(x)
                if ~isPlaying, break; end  % Stop if paused
                slider.Value = k;
                updatePlot(slider, []);
                drawnow;  % Force update
%                 pause(1/length(x));  % Adjust speed
            end
            isPlaying = false;
            playButton.String = 'Play';
        end
    end

    function for_every_plot(s)
        hold on;
        xlim([-x_w x_w]);
        ylim([-y_w y_w]);

        title(options.title);
        xlabel(["X coordinate (" options.units ")"])
        ylabel(["Y coordinate (" options.units ")"])

        if isnan(s) && ~strcmp(options.nan_message, "")
            text(options.ref(1)-100, options.ref(2), options.nan_message, 'FontSize', 20)
        end

        if strcmp(options.draw,"rect")
            if strcmp(options.pos, "center")
                l = options.ref(1) - options.width/2;
                b = options.ref(2) - options.height/2;
            end
            rectangle('Position',[l b options.width options.height], 'EdgeColor', 'k', 'LineWidth', 2)
        end

    end

end

function opt = set_options(var)
        % opt primeiramente guarda os valores padrões de cada variável
        % O ponto 'ref' vai ser ou o centro ou o canto (inferior esquerdo),
        % a depender de 'pos', do retângulo de referência 

        % -> Adicionar a opção de desenhar elipses
        % Para 'past' e 'thr', se 100, não altero os valores do ponto
        % atual até 100 pra trás, e mais pra trás começo a aplicar um dos
        % modos
        opt = struct('units', 'm', 'play_step', 1, 'draw', NaN, 'width', "auto", 'height', "auto", 'ref', [0,0], ...
            'pos', "center", 'title', "Interactive plot", 'color', 'r', ...
            'past', "keep", 'thr', NaN, 'nan_message', "");
        optionNames = fieldnames(opt);
        
        % Se for passada alguma entrada adicional, é alterado o valor
        % correspondente em opt
        for pair = reshape(var, 2, [])
            inpName = pair{1};
            if any(strcmp(inpName, optionNames))
                opt.(inpName) = pair{2};
            else
                error('%s is not a recognized parameter name', inpName);
            end
        end
    end
function emd_interactive_plot(eye_movs, t_range, data)
% eye_movs: obtido de mix_emd
% t_range:  ínicio e fim, em ms, do intervalo a ser plotado

mov_colors = [[219, 162, 4];
              [242, 67, 36];
              [242, 67, 36]; 
              [33, 194, 92];
              [126, 9, 235]
             ];

mov_colors = mov_colors./255;

mov_style = {
    mov_colors(1,:), ':',  'Blink/noise', 1.5;
    mov_colors(2,:), '-',  'Fixation',    2.;
    mov_colors(3,:), '--', 'Pursuit',     1.5;
    mov_colors(4,:), '-',  'Saccade',     2.;
    mov_colors(5,:), ':',  'PSO',         2.;
};

s_factor = .8;
l_factor = .7;
mov_colors_1 = desaturate_colors(mov_colors, s_factor, l_factor);

dim_mov_style = {
    mov_colors_1(1,:), ':',  'Blink/noise', 1.;
    mov_colors_1(2,:), '-',  'Fixation',    1.5;
    mov_colors_1(3,:), '--', 'Pursuit',     1.;
    mov_colors_1(4,:), '-',  'Saccade',     1.5;
    mov_colors_1(5,:), ':',  'PSO',         1.5;
};

l = length(mov_style);
f_ratio = eye_movs.exp_setup.frequency/1000;

if nargin < 3
    data = eye_movs.data.filt_lin;
end
lims = eye_movs.all_lims;
lims(3,:) = lims(3,:) + 1;

x = data(1,:);
y = data(2,:);
t_total = length(data)* f_ratio;
clear data;
if size(t_range, 1) == 1
    t_range = t_range.';
end

L = length(lims);

my_screen_size = get(0, 'ScreenSize');  % Get screen dimensions [1x4]
fig_width = my_screen_size(3) * 0.8;    % 80% of screen width
fig_height = my_screen_size(4) * 0.6;   % 60% of screen height
left = (my_screen_size(3) - fig_width) / 2;
bottom = (my_screen_size(4) - fig_height) / 2;
fig = figure('name', ['Eye Movement Detection: ', eye_movs.emd_method, ' Method' ], 'Position', [left, bottom, fig_width, fig_height]);
% Legenda para pontos de início e fim da trajetória
subplot(2, 1, 1); hold on; 
traj_legend = gobjects(2, 1);
traj_legend(1) = plot(NaN, NaN, 'k.', 'MarkerSize', 15, 'DisplayName', 'Start'); 
traj_legend(2) = line('XData', NaN, 'YData', NaN, ...
    'Marker', 's', ...          % 's' for square
    'MarkerSize', 7, ...       % Size of the square
    'MarkerFaceColor', 'm', ... % Magenta fill
    'MarkerEdgeColor', 'm', ... % Magenta border (optional)
    'LineStyle', 'none', ...    % Remove connecting lines
    'DisplayName', 'End');
legend(traj_legend, 'Location', 'northwest');

subplot(2, 1, 2); hold on;

% Create "dummy" lines for the legend (invisible data)
ev_types = unique(lims(3,:));
h_legend = gobjects(length(ev_types), 1);  % Handle array for legend entries
count = 1;
for k = 1:l
    if ismember(k, ev_types)
        h_legend(count) = plot(NaN, NaN, ...
            'Color', mov_style{k, 1}, ...
            'LineWidth', mov_style{k, 4}, ...
            'LineStyle', mov_style{k, 2}, ...
            'DisplayName', mov_style{k, 3});
        count = count + 1;
    end
end
legend(h_legend, 'Location', 'best');  % Add legend to current subplot

x_traces     = gobjects(1, L);
y_traces     = gobjects(1, L);
dim_x_traces = gobjects(1, L);
dim_y_traces = gobjects(1, L);
xy_traj      = gobjects(1, L + 2);
for i=1:L
    t_idx = lims(1,i):lims(2,i);    % Índices de tempo
    m_idx = lims(3,i);              % Índices do tipo de movimento
    subplot(2,1,1); hold on; 
    xy_traj(i) = plot(x(t_idx), y(t_idx), ...
                                  'Color', mov_style{m_idx, 1}, ...
                                  'LineWidth', mov_style{m_idx, 4}, ...
                                  'LineStyle', mov_style{m_idx, 2}, ...
                                  'HandleVisibility', 'off');
    
    subplot(2,1,2); hold on; 
    dim_x_traces(i) = plot((t_idx-1)*f_ratio, x(t_idx), 'Color', dim_mov_style{m_idx, 1}, ...
                       'LineWidth', dim_mov_style{m_idx, 4}, ...
                       'LineStyle', dim_mov_style{m_idx, 2}, ...
                       'HandleVisibility', 'off',...
                       'Visible', 'off');
    dim_y_traces(i) = plot((t_idx-1)*f_ratio, y(t_idx), 'Color', dim_mov_style{m_idx, 1}, ...
                       'LineWidth', dim_mov_style{m_idx, 4}, ...
                       'LineStyle', dim_mov_style{m_idx, 2}, ...
                       'HandleVisibility', 'off', ...
                       'Visible', 'off');
    x_traces(i) = plot((t_idx-1)*f_ratio, x(t_idx), 'Color', mov_style{m_idx, 1}, ...
                       'LineWidth', mov_style{m_idx, 4}, ...
                       'LineStyle', mov_style{m_idx, 2}, ...
                       'HandleVisibility', 'off',...
                       'PickableParts', 'visible',...
                       'ButtonDownFcn', @(src,evt) highlight_curve(fig, 'x'));
    y_traces(i) = plot((t_idx-1)*f_ratio, y(t_idx), 'Color', mov_style{m_idx, 1}, ...
                       'LineWidth', mov_style{m_idx, 4}, ...
                       'LineStyle', mov_style{m_idx, 2}, ...
                       'HandleVisibility', 'off', ...
                       'PickableParts', 'visible',...
                       'ButtonDownFcn', @(src,evt) highlight_curve(fig, 'y'));
end

traj_ax  = subplot(2,1,1); hold on; title('Eye position (x,y)'); xlabel(['x coordinate (' eye_movs.data.unit ')']); ylabel(['y coordinate (' eye_movs.data.unit ')'])
% Cria bola do inicio e quadrado do fim
xy_traj(L+1) = plot(x(1), y(1), 'k.', 'MarkerSize', 15, 'HandleVisibility', 'off'); 
xy_traj(L+2) = line('XData', NaN, 'YData', NaN, ...
    'Marker', 's', ...          % 's' for square
    'MarkerSize', 7, ...       % Size of the square
    'MarkerFaceColor', 'm', ... % Magenta fill
    'MarkerEdgeColor', 'm', ... % Magenta border (optional)
    'LineStyle', 'none', ...    % Remove connecting lines
    'DisplayName', 'End', 'HandleVisibility', 'off');

% Delimita o gráfico e a região de tela
screen_size = pixel_to_dva(eye_movs.exp_setup.resolution,"dist", 45,"pixel_dva_ratio", 27.854);
w = screen_size(1); h = screen_size(2);
rectangle('Position', [-w/2, -h/2, w, h], ...
          'EdgeColor', 'k', ...  % Red contour
          'LineWidth', 2, ...    % Thickness
          'LineStyle', '-');     % Solid line
ax_lims = [-eye_movs.data.pos_lims eye_movs.data.pos_lims];
set(traj_ax, 'XLim', ax_lims, 'XLimMode', 'manual', 'YLim', [ax_lims(1) + 5, ax_lims(2) - 5], 'YLimMode', 'manual', 'DataAspectRatio', [1 1 1]);

coord_ax = subplot(2,1,2); hold on; title('Eye coordinates over time'); xlabel(['Time (' eye_movs.time_unit ')']); ylabel(['Eye coordinates (' eye_movs.data.unit ')'])


%% Parte interativa
% Define o struct para ser repassado entre as funções
aux = struct();

aux.ax      = coord_ax;     % No eixo principal, ficam salvas as coordenadas
aux.ax2     = traj_ax;      % No eixo secundário, a trajetória
aux.h_xy    = xy_traj;

aux.fig     = fig;
aux.time    = ((1:t_total) - 1)*f_ratio;
aux.total_time = t_total;
aux.all_lims   = lims; aux.all_lims(3,:) = 1:L;

aux.x       = x;
aux.h_x     = x_traces;
aux.h_x_dim = dim_x_traces;

aux.y       = y;
aux.h_y     = y_traces;
aux.h_y_dim = dim_y_traces;

aux.current_start = t_range(1);
aux.winlength = t_range(2)-t_range(1);
aux.zoom_factor = 1.2;
aux.highlighted = 'none';
aux.ylim_cb_state  = 0;
aux.track_cb_state = 1;
aux.xlim = NaN;

clear eye_movs coord_ax traj_ax x x_traces dim_x_traces y y_traces dim_y_traces;


% Posição do ax para redimensioná-lo e posicionar os demais elementos
original_pos = get(aux.ax, 'Position');
pos_change = .07;
% Redefine a posição de ax
new_ax_pos = [original_pos(1), original_pos(2)+pos_change, original_pos(3), original_pos(4)-pos_change];

set(aux.ax, 'Position', new_ax_pos);

%% Define os objetos interativos
% (1) Botão de limpar na tela toda
set(fig, 'WindowButtonDownFcn', @(src,evt) clear_highlights(fig));

% (2) Slider
slider_pos = [original_pos(1), original_pos(2)-pos_change, original_pos(3), .5*pos_change];
aux.slider = uicontrol('Style', 'slider',...
    'Parent', fig,...
    'Min', 0, 'Max', t_total-(t_range(2) - t_range(1)),...   % O máximo é o total menos o tamanho da janela
    'Value', aux.current_start,...
    'Units','normalized',...
    'Position', slider_pos,...
    'Callback', @(src,evt) update_plot(fig)...
                           );
%                            'Callback', @(src,evt) update_plot(fig)

% (3) Botões de zoom
button_size = pos_change/2;
button_dist = .01;
plus_button_pos  = [slider_pos(1)+slider_pos(3)+button_dist, slider_pos(2)+.5*(button_size+button_dist), button_size, button_size];
minus_button_pos = [slider_pos(1)+slider_pos(3)+button_dist, slider_pos(2)-.5*(button_size+button_dist), button_size, button_size];

uicontrol('Style', 'pushbutton',...
    'Parent', fig,...
    'String', '-',...
    'Units','normalized',...
    'Position', minus_button_pos,...
    'Callback', @(src,evt) zoom_out(fig));
    
uicontrol('Style', 'pushbutton',...
    'Parent', fig,...
    'String', '+',...
    'Units','normalized',...
    'Position', plus_button_pos,...
    'Callback', @(src,evt) zoom_in(fig));
    

aux.top_button_pos = [(original_pos(1)+original_pos(3)+button_dist) (original_pos(2)+.5*(original_pos(4)+button_dist)) button_size button_size];
aux.bot_button_pos = [(original_pos(1)+original_pos(3)+button_dist) (original_pos(2)-button_size+.5*(original_pos(4)-button_dist)) button_size button_size];

% (4) Botões de destaque (highlight)
aux.btn_x = uicontrol('Style', 'pushbutton',...
      'Parent', fig,...
      'String', 'x(t)',...
      'Units','normalized',...
      'Position', aux.top_button_pos,...
      'Callback', @(src,evt) highlight_curve(fig, 'x'));

aux.btn_y = uicontrol('Style', 'pushbutton',...
      'Parent', fig,...
      'String', 'y(t)',...
      'Units','normalized',...
      'Position', aux.bot_button_pos,...
      'Callback', @(src,evt) highlight_curve(fig, 'y'));

vrt_ax_button_pos = [slider_pos(1), slider_pos(2)-.5*button_dist-button_size, .25*slider_pos(3), button_size];
% (5) Botão para definir ylim
aux.btn_vert_axis = uicontrol('Style', 'checkbox',...
    'Parent',fig,...
    'String','Auto-adjust y axis',...
    'Units','normalized',...
    'Position', vrt_ax_button_pos,...
    'Callback', @(src,~) auto_adjust_y(src, fig));

track_button_pos = [slider_pos(1) + .5*slider_pos(3), slider_pos(2)-.5*button_dist-button_size, .25*slider_pos(3), button_size];
% (6) Botão para acompanhar os limites 
aux.btn_keep_track = uicontrol('Style', 'checkbox',...
    'Parent',fig,...
    'String','Track start and end',...
    'Units','normalized',...
    'Position', track_button_pos,...
    'Callback', @(src,~) track_traj_lims(src, fig));

% Pode parecer confuso, mas a ideia é que plote início e fim corretamente
% usando aux.track_cb_state = 1 mas que depois...
guidata(fig, aux); update_plot(fig);
% ... Seja definido 0 como valor inicial (sem tracking)
aux.track_cb_state = 0; guidata(fig, aux);

end

%% Funções de GUI

% Atualização do gráfico
function update_plot(fig)
    aux = guidata(fig);
    if strcmp(get(aux.slider, 'Enable'), 'on')
        aux.current_start = get(aux.slider, 'Value');
    else
        aux.current_start = 0;  % Default when slider is disabled
    end
    xlim(aux.ax, [aux.current_start, aux.current_start+aux.winlength]);
    aux.xlim = [aux.current_start, aux.current_start+aux.winlength];

    
    % Toggle ylim behavior
    if aux.ylim_cb_state  == 1
        % Allow auto-adjustment when checked
        set(aux.ax, 'YLimMode', 'auto');
    else
        % Keep current ylim when unchecked
        ylim(aux.ax, get(aux.ax, "ylim"));
    end

    
    
    guidata(fig, aux);
    update_button_positions(fig);
    if aux.track_cb_state == 1
        adjust_trajectory(fig);
    end
end

% Atualiza a posição dos botões 
function update_button_positions(fig)
    aux = guidata(fig);
    x_right_edge = aux.current_start + aux.winlength;
    [~, idx] = min(abs(aux.time - x_right_edge));
    
    h_x = aux.x(idx);
    h_y = aux.y(idx);
    if h_x >= h_y
        set(aux.btn_x, 'Position', aux.top_button_pos);
        set(aux.btn_y, 'Position', aux.bot_button_pos);
    else
        set(aux.btn_x, 'Position', aux.bot_button_pos);
        set(aux.btn_y, 'Position', aux.top_button_pos);
    end
end

% Funções de zoom
function zoom_in(fig)
    aux = guidata(fig);
    aux.winlength = max(aux.winlength/aux.zoom_factor, 0.1);
    set(aux.slider, 'Enable', 'on');
    set(aux.slider, 'Max', aux.total_time-aux.winlength);
    if get(aux.slider, 'Value') > get(aux.slider, 'Max')
        set(aux.slider, 'Value', get(aux.slider, 'Max'));
    end
    guidata(fig, aux);
    update_plot(fig);
end

function zoom_out(fig)
    aux = guidata(fig);
    aux.winlength = min(aux.winlength*aux.zoom_factor, aux.total_time);

    % Check if fully zoomed out
    if aux.winlength >= aux.total_time
        aux.winlength = aux.total_time;  % Ensure exact match
        set(aux.slider, 'Enable', 'off');  % Disable slider (no effect when fully zoomed)
        set(aux.slider, 'Value', 0);       % Reset to avoid warning (value irrelevant)
    else
        set(aux.slider, 'Enable', 'on');
        set(aux.slider, 'Max', aux.total_time - aux.winlength);
    end

    set(aux.slider, 'Value', min(get(aux.slider, 'Max'), get(aux.slider, 'Value')));
    guidata(fig, aux);
    update_plot(fig);
end

% Funções de destaque (highlight)
function highlight_curve(fig, curve)
    aux = guidata(fig);
    
    % If clicking the already highlighted curve, do nothing
    if strcmp(aux.highlighted, curve)
        return;
    end
    
    % Update highlight state
    aux.highlighted = curve;
    
    % Update visuals
    switch curve
        case 'x'
            set(aux.h_x, 'Visible', 'on');
            set(aux.h_y, 'Visible', 'off');
            set(aux.h_x_dim, 'Visible', 'off');
            set(aux.h_y_dim, 'Visible', 'on');
        case 'y'
            set(aux.h_x, 'Visible', 'off');
            set(aux.h_y, 'Visible', 'on');
            set(aux.h_x_dim, 'Visible', 'on');
            set(aux.h_y_dim, 'Visible', 'off');
    end
    
    guidata(fig, aux);
end

function clear_highlights(fig)
    aux = guidata(fig);
    aux.highlighted = 'none';
    set(aux.h_x, 'Visible', 'on');
    set(aux.h_y, 'Visible', 'on');
    set(aux.h_x_dim, 'Visible', 'off');
    set(aux.h_y_dim, 'Visible', 'off');
    guidata(fig, aux);
end

function auto_adjust_y(src, fig)
    % Get current axes and checkbox state
    aux = guidata(fig);
    is_checked = get(src, 'Value'); % 1=checked (auto), 0=unchecked (manual)

    aux.ylim_cb_state = is_checked;
    guidata(fig, aux);
    update_plot(fig);
end

function adjust_trajectory(fig)
    aux = guidata(fig);

    compl_xlim = trim_overlapping_events([0; aux.total_time], aux.xlim');
    lims_to_remove = trim_overlapping_events_v2(aux.all_lims, aux.xlim'); 
    if ~isempty(lims_to_remove)
        evs_to_remove = lims_to_remove(3,:); clear lims_to_remove;
    else
        evs_to_remove = [];
    end
    lims_to_keep = trim_overlapping_events_v2(aux.all_lims, compl_xlim);  evs_to_keep = lims_to_keep(3,:);
    start_time = aux.all_lims(1,lims_to_keep(3,1)); 
    end_time   = aux.all_lims(2, lims_to_keep(3, end));
    clear lims_to_keep;
    
    set(aux.h_xy(evs_to_remove), 'Visible', 'off');
    set(aux.h_xy(evs_to_keep), 'Visible', 'on');
    set(aux.h_xy(end-1), 'XData', aux.x(start_time), 'YData', aux.y(start_time));
    set(aux.h_xy(end), 'XData', aux.x(end_time), 'YData', aux.y(end_time));
end

function track_traj_lims(src, fig)
    % Get current axes and checkbox state
    aux = guidata(fig);
    is_checked = get(src, 'Value'); % 1=checked (track), 0=unchecked (not track)

    aux.track_cb_state  = is_checked;
    guidata(fig, aux);
    update_plot(fig);
end
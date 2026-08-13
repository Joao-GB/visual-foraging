function scrollable_plot(plot_handle, total_time, plot_window)
    % ADD_SCROLLABILITY - Adds scrolling functionality to existing figure/axes
    % plot_handle:   handle to figure or axes
    % total_time:    total time of the curve plottted
    % plot_window:   initial window to display (in seconds)


    % Determine input type and get relevant axes
    if strcmp(get(plot_handle, 'Type'), 'figure')
        % If input is a figure, find all axes
        all_axes = findobj(plot_handle, 'Type', 'axes');
        
        % Exclude legends and colorbars
        valid_axes = [];
        for ax = all_axes'
            if isempty(findobj(ax, 'Type', 'legend')) && isempty(findobj(ax, 'Tag', 'Colorbar'))
                valid_axes = [valid_axes; ax];
            end
        end
        
        % Use the first axes if only one exists, otherwise use the second
        if length(valid_axes) >= 2
            ax = valid_axes(2); % Target subplot(2,1,2)
        else
            ax = valid_axes(1); % Fallback to the only axes
        end
    else
        % Input is already an axes handle
        ax = plot_handle;
    end
    
    % Get parent figure
    fig = ancestor(ax, 'figure');
    
    % Store original axes position for proper UI element placement
    original_pos = get(ax, 'Position');
    
    % Adjust axes position to make room for controls
    pos_change = .07;
    new_ax_pos = [original_pos(1), original_pos(2)+pos_change, original_pos(3), original_pos(4)-pos_change];
    set(ax, 'Position', new_ax_pos);
    
    % Initialize parameters
    current_start = plot_window(1);
    
    % Create slider
    slider_pos = [original_pos(1), original_pos(2)-pos_change, original_pos(3), .5*pos_change];
    uicontrol('Style', 'slider',...
        'Parent', fig,...
        'Min', 0, 'Max', total_time-(plot_window(2) - plot_window(1)),...   % O máximo é o total menos o tamanho da janela
        'Value', current_start,...
        'Units','normalized',...
        'Position', slider_pos,...
        'Callback', @update_plot);
    
    button_size = pos_change/2;
    button_dist = .01;
    % Create zoom buttons
    uicontrol('Style', 'pushbutton',...
        'Parent', fig,...
        'String', '-',...
        'Units','normalized',...
        'Position', [slider_pos(1)-(button_dist + button_size), slider_pos(2), button_size, button_size],...
        'Callback', @zoom_out);
    
    uicontrol('Style', 'pushbutton',...
        'Parent', fig,...
        'String', '+',...
        'Units','normalized',...
        'Position', [slider_pos(1)+slider_pos(3)+button_dist, slider_pos(2), button_size, button_size],...
        'Callback', @zoom_in);
    
    % Store data in axes userdata for callbacks

    set(ax, 'UserData', struct(...
        'winlength', (plot_window(2)-plot_window(1)),...
        'total_time', total_time, ...
        'zoom_factor', 1.2));
    
    % Set initial view
    update_plot();
    
    % Callback functions
    function update_plot(~,~)
        ud = get(ax, 'UserData');
        slider = findobj(fig, 'Style', 'slider');

        if strcmp(get(slider, 'Enable'), 'on')
            current_start = get(slider, 'Value');
        else
            current_start = 0;  % Default when slider is disabled
        end

        current_ylim = get(ax, "ylim");
%         if current_start+ud.winlength > ud.total_time
%             current_start = ud.total_time - ud.winlength;
%         end
        xlim(ax, [current_start current_start+ud.winlength]);
        ylim(ax, current_ylim);
    end

    function zoom_in(~,~)
        ud = get(ax, 'UserData');
        ud.winlength = max(ud.winlength/ud.zoom_factor, 0.1);
        slider = findobj(fig, 'Style', 'slider');
        set(slider, 'Enable', 'on');
        set(slider, 'Max', ud.total_time-ud.winlength);
        if get(slider, 'Value') > get(slider, 'Max')
            set(slider, 'Value', get(slider, 'Max'));
        end
        set(ax, 'UserData', ud);
        update_plot();
    end

    function zoom_out(~,~)
        ud = get(ax, 'UserData');
        ud.winlength = min(ud.winlength*ud.zoom_factor, ud.total_time);

        slider = findobj(fig, 'Style', 'slider');

        % Check if fully zoomed out
        if ud.winlength >= ud.total_time
            ud.winlength = ud.total_time;  % Ensure exact match
            set(slider, 'Enable', 'off');  % Disable slider (no effect when fully zoomed)
            set(slider, 'Value', 0);       % Reset to avoid warning (value irrelevant)
        else
            set(slider, 'Enable', 'on');
            set(slider, 'Max', ud.total_time - ud.winlength);
        end

        slider_max = get(slider, 'Max'); slider_value = get(slider, 'Value');
        set(slider, 'Value', min(slider_max, slider_value));
        set(ax, 'UserData', ud);
        update_plot();
    end
end
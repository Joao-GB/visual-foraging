function d = eye_filter(data, fs, window, fc, order, type)
% Aplica filtro de mediana e Butterworth para suavizar os dados de
% eye-tracking (mas não remove grandes amplitudes, por isso usar depois da
% remoção de blinks)

    % fs    : frequência de amostragem dos dados
    % window: janela usada para medfilt1
    % fc    : cutoff frequency a ser passada para o filtro
    % order : ordem do filtro Butterworth
    % type  : tipo do filtro Butterworth (e.g., 'high' ou 'low')

    % Aplica um filtro de mediana com janela +-25 ms -- como seus dados
    % estão a 120Hz, usa janela de tamanho +-3. No nosso caso, a janela deve
    % ser de tamanho +-25
    h_win = ceil(window/2);

    if window > 0

        x_pad = apply_padding(double(data(1,:)), h_win, 'add', 'symmetric');
        y_pad = apply_padding(double(data(2,:)), h_win, 'add', 'symmetric');
        x = apply_padding(medfilt1(x_pad, window), h_win, 'remove');
        y = apply_padding(medfilt1(y_pad, window), h_win, 'remove');
    else
        x = double(data(1,:));
        y = double(data(2,:));
    end


    % Filtro passa-alto para remover ruído de baixa frequência
%     [b_hp, a_hp] = butter(order, 0.5/(fs/2), 'high');
%     x = filtfilt(b_hp, a_hp, x);
%     y = filtfilt(b_hp, a_hp, y);


    [b, a] = butter(order,fc/(fs/2),type); % Butterworth filter
                                           % retorna os coeficientes da transfer function do filtro
    % Aplicação com ida e volta (filtfilt) para fazê-lo não-causal
    d(1, :) = filtfilt(b, a, x);
    d(2, :) = filtfilt(b, a, y);

end
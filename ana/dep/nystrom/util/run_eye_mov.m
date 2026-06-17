switch_project('mestrado')
d541 = load("C:\Users\joaog\OneDrive\Área de Trabalho\USP\Mestrado\Projeto\data\free view\541\Anno_541_export.mat");
pre_data = d541.data.EYES; fs = 1000;

% Primeira tentativa
st = eye_mov_detection(pre_data, fs, 'pixel');

% Usando detecção de Ito
st1 = ito_grun_emd(pre_data, fs, 'pixel');

post_data(1:2,:) = st.saccades.data;
post_data(3:4,:) = st.fixations.data;
post_data(5:6,:) = st.pursuits.data;

% Plotar apenas sacadas
plot_at_eeglab(post_data(1:2, :))

plot(post_data(1:2, 6000:14000).')
title("Candidate Saccades")
xlabel("Time (ms)")
ylabel("Eye position (dva)")

%Plotar apenas fixações
plot_at_eeglab(post_data(3:4, :));
% Plotar apenas pursuits
plot_at_eeglab(post_data)

% Obter largura da tela em graus
len = pixel_to_dva([1680, 1050], 'dist', 45', 'pixel_dva_ratio', 27.854);

% Para plotar a trajetória interativa
d = st.nan_data;
interactive_trajectory(d, fs, 'units', 'dva', 'title', "Eye trace", 'play_step', 40, 'draw', "rect", 'width', len(1), 'height', len(2), 'past', "lighter", 'thr', 2000, 'nan_message', "blink/noise");
function fakeLoadingScreen(tkP, dpP, drP, prm, loadingMode, txP, ori)
    if nargin < 5, loadingMode = 'opening'; txP = []; ori = []; end
    % Tela de 'Carregando' falsa com barra e dicas
    aux = startFake(tkP, dpP, drP, prm, loadingMode, txP, ori);
    endFake(aux, dpP, drP, prm)
end
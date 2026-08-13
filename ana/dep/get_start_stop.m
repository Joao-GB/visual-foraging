function [start, stop] = get_start_stop(bol)
% Encontra em um booleano onde começa e onde termina cada intervalo de 1s
 
    % Se for vetor de índices em vez de booleano, é convertido
    if ~islogical(bol)
        aux = false(1, max(bol)); aux(bol) = true;
        bol = aux; clear aux;
    end
    l = length(bol);
    
    d_bol = diff(bol);      % Retorna aux(2:end) - aux(1:end-1)
    start  = find(d_bol == 1) + 1; 
    stop = find(d_bol == -1);
    
    % Se, por algum azar, os registros terminarem (ou começarem) durante o uma
    % evento de interesse, adiciono fim (ou início) do último (primeiro) evento
    % -> Se continuarem com tamanhos diferentes, há algo de errado com os dados...
    if bol(1) == 1
        start = [1 start];
    end
    if bol(end) == 1
        stop = [stop l];
    end
end
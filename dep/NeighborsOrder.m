function [nbhd] = NeighborsOrder(points, Pidx)
    % points: matriz Nx2 dos pontos de interesse
    % Pidx:   indice de um dos pontos de points cuja vizinhança se quer
    %          ordenar

    if size(points,2) > size(points,1), points = points'; end
    
    % 1. Obtém as coordenadas de P
    P = points(Pidx,:);
    
    % 2. Calcula a distancia entre P e os demais
    qDists = sum((points - P).^2, 2);
    [~, nbhd] = sort(qDists);
    nbhd = nbhd';

    % Remove o próprio ponto
    nbhd(nbhd == Pidx) = [];
end
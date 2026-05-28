function [nbhd] = NeighborsOrder(points, Pidx, P1idx)
    % points: matriz Nx2 dos pontos de interesse
    % Pidx:   indice de um dos pontos de points cuja vizinhança se quer
    %          ordenar
    % P1idx:  índice do ponto alvo da sacada
    if size(points,2) > size(points,1), points = points'; end
    if isempty(P1idx), P1idx = Pidx; end
    
    % 1. Obtém as coordenadas de P
    P  = points(Pidx,:);
    
    % 2. Calcula a distancia real entre P e os demais (adicionado sqrt)
    dists = sqrt(sum((points - P).^2, 2));
    
    % 3. Calcula a diferença absoluta para a distância alvo (P até P1)
    targetDist = dists(P1idx);
    distDiff = abs(dists - targetDist);
    
    % 4. Ordena pela menor diferença (mais parecida)
    [~, nbhd] = sort(distDiff);
    nbhd = nbhd';
    
    % Remove os próprios pontos
    nbhd(nbhd == Pidx) = [];
    nbhd(nbhd == P1idx) = [];
end
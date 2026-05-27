function [nbhd1, nbhd2, nbhdElse] = getHexNeighborhoods2(points, Pidx, nbDist)
    if size(points,2) > size(points,1), points = points'; end

    S = size(points,1);
    if isempty(Pidx) || Pidx <= 0 || Pidx > S
        aux = 1:S;
        nbhd1 = aux; nbhd2 = aux; nbhdElse = aux;
        return
    end
    
    % points: matriz Nx2 dos pontos de interesse
    % Pidx:   indice de um dos pontos de points cuja izinhança se quer particionar
    
    % 1. Obtém as coordenadas de P
    P = points(Pidx,:);
    
    % 2. Calcula a distancia entre P e os demais
    dists = sqrt(sum((points - P).^2, 2));
    
    % Sanity check: garante ao menos 2 pontos vizinhos
    if length(dists) < 2
        nbhd1 = []; nbhd2 = []; nbhdElse = [];
        return;
    end
    
    % 4. Classifica vizinhança conforme a banda de distancia
    
    % Limiares entre ambas as vizinhanças possíveis
    limit1 = 1.4 * nbDist;
    limit2 = 2.5 * nbDist;
    
    % Nbhd1: vizinhos de primeira ordem (sem contar o próprio P)
    nbhd1 = find(dists > eps & dists < limit1)';
    
    % Nbhd2: vizinhos de segunda ordem
    nbhd2 = find(dists > limit1 & dists < limit2)';
    
    % NbhdElse: todos os demais vizinhos, distantes
    nbhdElse = find(dists > limit2)';
end
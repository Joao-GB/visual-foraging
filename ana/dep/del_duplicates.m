% Substitui todas subsequências de S que tenham elementos iguais por um
% único elemento, deixando-a sem repetições consecutivas
%
% INPUT
% S:        matriz (array) que contenha a sequência de símbolos
%
% OUTPUT
% T:        sequência sem repetições consecutivas

% Date:   10/2023
% Author: João Borges


function T = del_duplicates(S)
% Cria uma matriz que vale 1 onde há repetição ao comparar termos
% consecutivos. Nunca exclui o primeiro termo, pois se for igual ao
% próximo, este que é removido.
aux = [logical(0) S(2:end) == S(1:end-1)];
T = S(~aux);
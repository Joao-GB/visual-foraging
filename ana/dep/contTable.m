function [T, I, d, c] = contTable(E, O, showTable)
    if nargin < 3, showTable = 0; end 
    T = zeros(2,2); I = cell(2,2);
    
    % Linhas (1: presente, 2: ausente) | Colunas (1: sim, 2: não)
    
    % Hit (presente + sim)
    I{1,1} = (E == 1 & O == 1);
    T(1,1) = sum(I{1,1});
    
    % Miss (presente + não)
    I{1,2} = (E == 1 & O == 0);
    T(1,2) = sum(I{1,2});
    
    % Alarme falso (ausente + sim)
    I{2,1} = (E == 0 & O == 1);
    T(2,1) = sum(I{2,1});
    
    % Rejeição correta (ausente + não)
    I{2,2} = (E == 0 & O == 0);
    T(2,2) = sum(I{2,2});
    
    % Cálculo das taxas
    nPresent = T(1,1) + T(1,2); % Sianis 
    nAbsent  = T(2,1) + T(2,2); % Total Signal Absent trials
    
    H_rate  = T(1,1) / nPresent;
    FA_rate = T(2,1) / nAbsent;
    
    % Correção de Hautus para evitar infinitos nas operações
    if H_rate == 1,   H_rate = 1 - (0.5 / nPresent); end
    if H_rate == 0,   H_rate = 0.5 / nPresent;       end
    if FA_rate == 1, FA_rate = 1 - (0.5 / nAbsent);  end
    if FA_rate == 0, FA_rate = 0.5 / nAbsent;        end
    
    % Métricas de SDT
    zH  = norminv(H_rate);
    zFA = norminv(FA_rate);
    d = zH - zFA;
    c = -(zH + zFA) / 2;
    
    if showTable
        fprintf('\n======================================\n');
        fprintf('             RESPOSTA DO SUJEITO\n');
        fprintf('              Sim (1)       Não (0)\n');
        fprintf('            -----------   -----------\n');
        fprintf('Alvo  (1)   %4d (H)      %4d (M)\n',  T(1,1), T(1,2));
        fprintf('Distr (0)   %4d (FA)     %4d (CR)\n', T(2,1), T(2,2));
        fprintf('--------------------------------------\n');
        fprintf('Totais: Alvo = %d, Distrator = %d\n', nPresent, nAbsent);
        fprintf('Taxas:   Hit = %.2f%%,    FA = %.2f%%\n', H_rate*100, FA_rate*100);
        fprintf('          d'' = %.4f,     c = %.4f\n', d, c);
    end
end
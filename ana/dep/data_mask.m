function data = data_mask(data, mask, m_type)
% Assume que os dados são 2D, com cada linha representando uma atividade ao
% longo do tempo e que a máscara será aplicada sobre os mesmos
% subintervalos para todas as atividade

% Posso acrescentar outros tipos de máscara, como pchip
    data(:, mask) = NaN;
    [L,~] = size(data);
    if exist('m_type', 'var')
        for i=1:L; data(i,:) = fillmissing(data(i,:), m_type); end
    end
end
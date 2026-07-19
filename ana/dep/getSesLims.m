function sesLims = getSesLims(mat, messages, expIdx)
%% OBS: entenda sessões dentro do arquivo como subsessões do par (sujeito,sessão)
% Pode haver mais de um início de (sub)sessão experimental por arquivo, mas
% no máximo uma completa
    sesLims = [];

    expSesOnMsg = sprintf(mat.prm.msg.on.ses{1}, mat.prm.msg.suffix{expIdx});
    expSesOnIdx = [find(contains(messages, expSesOnMsg))];

   if isempty(expSesOnIdx), return; end

    auxIdx = [expSesOnIdx length(messages)];
    expSesOffIdx = zeros(1, length(expSesOnIdx));
    for i=1:length(expSesOnIdx)
        % Procura, entre cada par de mensagem de início de sessão,
        % uma mensagem de fim de sessão
        offAuxIdx = [find(contains(messages(auxIdx(i):auxIdx(i+1)), mat.prm.msg.off.ses{3}))] + expSesOnIdx(i) - 1;
        if length(offAuxIdx) ~= 1
            warning('Mais de um fim associado ao início da sessão  %d do arquivo! Tomando o primeiro', i)
            expSesOffIdx(i) = offAuxIdx(1);
        else
            expSesOffIdx(i) = offAuxIdx;
        end
    end
    % Exclui as sessões 
    expSesOnIdx(expSesOffIdx == 0)  = [];
    expSesOffIdx(expSesOffIdx == 0) = [];

    %% Encontra a última sessão experimental completa
    sesCpMsg = sprintf(mat.prm.msg.off.ses{1}, mat.prm.msg.suffix{expIdx});
    auxIdx = [find(strcmp(messages(expSesOffIdx), sesCpMsg), 1, "last")];
    if isempty(auxIdx)
        warning('Não há sessões experimentais concluídas, usando a última interrompida');
        sesCpMsg = sprintf(mat.prm.msg.off.ses{2}, mat.prm.msg.suffix{3});
        auxIdx = min([length(messages), [find(strcmp(messages(expSesOffIdx), sesCpMsg), 1, "last")]]);
    end

    sesOnIdx  = expSesOnIdx(auxIdx);
    sesOffIdx = expSesOffIdx(auxIdx);
    sesLims   = [sesOnIdx; sesOffIdx];
    clear expSesOnMsg offAuxIdx expSesOnIdx expSesOffIdx auxIdx i;
end
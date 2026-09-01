function outTrl = foragingPSApermute(trl)
% Apenas troca tanto categorias, respostas e hits entre condições sacádica
% e não sacádica (por mais que getPSAeffect não use os hits)
% Precisa ver de trocar mais coisa dependendo das funções que forem usar


    K = numel(trl);
    outTrl = trl;

    swap = rand(1,K) < 0.5;

    for k = 1:K
        if swap(k)
            outTrl(k).probeCat       = trl(k).nSaccProbeCat;
            outTrl(k).nSaccProbeCat  = trl(k).probeCat;
            outTrl(k).probeResp      = trl(k).nSaccProbeResp;
            outTrl(k).nSaccProbeResp = trl(k).probeResp;
        end
    end
end
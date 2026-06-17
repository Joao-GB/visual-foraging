function b = nystrom_get_blinks(data, pupil, vel, acc, options)
% NYSTROM_GET_BLINKS identifica piscadas e ruído com base no diâmetro da
% pupila, coordenadas oculares e velocidade e aceleração dos olhos,
% conforme descrito em Nystrom & Holmqvist (2010). 
% INPUT:
% data   : posição dos olhos em função do tempo, em dva (cf. pixel_to_dva)
% pupil  : tamanho da pupila em função do tempo
% options: struct com as seguintes opções
%   posThr: limiar de posição, em dva, a partir do qual se considera o 
%           movimento não fisiológico ou se descarta por sair de uma área
%           de interesse
%   pupThr: limiar do tamanho da pupila, em ??, abaixo do qual se considera
%           que os olhos estão fechados
%   velThr: limiar de velocidade, em dva/s, a partir do qual se considera
%           o movimento não fisiológico
%   accThr: limiar de aceleração, em dva/s^2
%   itvThr: intervalo mínimo entre duas piscadas, em pts, para que sejam
%           mantidas separadas. Se o intervalo entre piscadas consecutivas 
%           for muito curto, fundo-as, pois não deve contribuir para a 
%           compreensão do comportamento ocular
%
% OUTPUT:
% b       : struct com informações sobre os intervalos de piscada/ruído

    arguments
        data  {mustBeNumeric} = [];
        pupil {mustBeNumeric} = [];
        vel   {mustBeNumeric} = [];
        acc   {mustBeNumeric} = [];
        options.posThr  {mustBeNumeric} = 50 
        options.velThr  {mustBeNumeric} = 1000;
        options.accThr  {mustBeNumeric} = 100000;
        options.pupThr                  = -1;
        options.itvThr  {mustBeNumeric} = 10;
    end
    
    pos_thr = options.posThr; pup_thr = options.pupThr;
    vel_thr = options.velThr; acc_thr = options.accThr;
    min_itv = options.itvThr;
    
    x = data(1, :); y = data(2, :);

% 1) Encontra os índices em que ao menos uma das coordenadas excede seu
    % respectivo limiar
    % (a) Limiar de posição
    if length(pos_thr) == 1
        aux = abs(sqrt(x.^2+y.^2)) > pos_thr;
    else
        aux = x > pos_thr(1) | y > pos_thr(2);
    end

    % (b) Limiares de velocidade e aceleração
    aux = aux | vel > vel_thr | acc > acc_thr;

    % (c) Limiar de tamanho de pupila
    if ~isempty(pupil)
        if pup_thr == -1; pup_thr = median(pupil)-3*std(pupil); end
        aux = aux | pupil <= max(pup_thr,0);
    end

% 2) Recalcula o início e o fim dos intervalos de piscadas
    [b_on, b_off] = get_start_stop(aux);

    % (a) Usa critério de velocidade mediana para decidir se cresce ou não
    %     a atual piscada
    med_vel = median(vel);
    L = length(x);
    for i = 1:length(b_on)
        while b_on(i)  - 1 >= 1 && vel(b_on(i)  - 1) >= med_vel; b_on(i)  = b_on(i)  - 1; end
        while b_off(i) + 1 <= L && vel(b_off(i) + 1) >= med_vel; b_off(i) = b_off(i) + 1; end
    end

    % (b) Se b_off de uma piscada anterior cruzar b_on da piscada seguinte,
    % são o mesmo evento
    b_aux = merge_overlapping_events([b_on; b_off]);

    % (c) Funde piscadas próximas
    b_aux = merge_adjacent_events(b_aux, min_itv);

    % (d) Obtém struct final
    b = mov_props(b_aux);
end
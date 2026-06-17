function [fLI, fLA, fTI, fD, fEP, fCV, detail] = getFixLimsOLD(e, et, times, pLims, fEv, fLims, eyePos, mode)
    % O modo indica se a fixação será procurada com relação ao fim ou ao
    % início do intervalo
    fLI=[]; fLA=[]; fTI=[]; fD=[]; fEP=[]; fCV = []; detail=[];
    %% 1 se quiser lidar com início de fixação; 0 para fins (P3, pois fixação começou na P2)
    if nargin < 8, mode = 1; end
    if mode
        if numel(pLims) == 1 % Procura a fixação mais próxima a começar antes
            aux = find(strcmp(e(1:pLims), fEv{1}), 1, "last");
        else                % Procura todas as fixações no intervalo, 
                aux = find(strcmp(e(pLims(1):pLims(2)), fEv{1})) + pLims(1) - 1;
        end
        N = numel(aux);
        if N == 0, return;
        else
            aux2 = find(strcmp(e(aux(1):end), fEv{2}),N) + aux(1) - 1;
            if numel(aux) ~= numel(aux2), return;
            else
                fLI = [aux; aux2];
            end
        end
    else
        if numel(pLims) == 1 % Procura a fixação mais próxima a terminar antes
            aux = find(strcmp(e(1:pLims), fEv{2}), 1, "last");
        else                % Procura todas as fixações no intervalo, 
            aux = find(strcmp(e(pLims(1):pLims(2)), fEv{2})) + pLims(1) - 1;
        end
        N = numel(aux);
        if N == 0, return;
        else
            aux2 = find(strcmp(e(1:aux(N)), fEv{1}),N, "last");
            if numel(aux2)~=numel(aux), return;
            else
                fLI = [aux2; aux];
            end
        end
    end


    fLA = reshape(et(:, fLI(2,:)), [2,N]);
    fTI = [arrayfun(@(x) find(times==x,1), fLA(1,:)); arrayfun(@(x) find(times==x,1), fLA(2,:))];
    fD = double(fLA(2,:)-fLA(1,:))/1000;

    fEP = zeros(2, size(fTI,2)); fCV = zeros(2, 2, size(fTI,2));
    for i=1:size(fTI,2)
        [fEP(:,i), fCV(:,:,i)] = fixStats(eyePos(:,fTI(1,i):fTI(2,i)));
    end

    % Agrega as fixações conforme os limites esperados (e.g., se houver
    % duas fixações dentro de um mesmo estímulo, considera como uma só

    nFix = size(fLims, 2);
    fLims(1,1) = 1; fLims(2,end) = numel(e);
    auxfLI = {}; auxfLA = {}; auxfTI = {}; auxfD = {}; auxfEP = {}; auxfCV = {};
    rfLI = zeros(2,nFix); rfLA = zeros(2,nFix); rfTI = zeros(2,nFix); rfD = zeros(1,nFix); rfEP = zeros(2,nFix); rfCV = zeros(2,2,nFix);

    for i=1:nFix
        % talvez não faça sentido incluir todas que comecem ou terminem
        idx = find((fLI(1,:) >= fLims(1,i) & fLI(1,:) <= fLims(2,i)) | (fLI(2,:) >= fLims(1,i) & fLI(2,:) <= fLims(2,i)));
        if isempty(idx)
            idx = find((fLI(1,:) <= fLims(1,i) & fLI(2,:) >= fLims(1,i)) | (fLI(1,:) <= fLims(2,i) & fLI(2,:) >= fLims(2,i)));
        end
        if isempty(idx), fLI = []; fLA = []; fTI = []; fD =[]; fEP = []; fCV = []; detail = []; return; end

        auxfLI = [auxfLI {fLI(:,idx)}]; %#ok<*AGROW> 
        auxfLA = [auxfLA {fLA(:,idx)}];
        auxfTI = [auxfTI {fTI(:,idx)}];
        auxfD = [auxfD {fD(idx)}];
        auxfEP = [auxfEP {fEP(:,idx)}];
        auxfCV = [auxfCV {fCV(:,:,idx)}];

        rfLI(:,i) = [fLI(1,idx(1));fLI(2,idx(end))];
        rfLA(:,i) = [fLA(1,idx(1));fLA(2,idx(end))];
        rfTI(:,i) = [fTI(1,idx(1));fTI(2,idx(end))];
        rfD(i) = double(rfLA(2,i)-rfLA(1,i))/1000;
        auxD = auxfD{end}./sum(auxfD{end}(:));
        rfEP(:,i) = auxfEP{end}*auxD';

        K = length(idx);
        n = auxfD{end}(:)';
        outfCV = zeros(2,2); outfEP = rfEP(:,i);
        
        for j = 1:K
            
            mu_i = auxfEP{end}(:,j);
            Sigma_i = auxfCV{end}(:,:,j);
            
            diff = mu_i - outfEP;
            
            outfCV = outfCV + ...
                (n(j)-1)*Sigma_i + ...
                n(j)*(diff*diff');
        end
        rfCV(:,:,i) = outfCV;
    end

    detail.fLI = auxfLI; detail.fLA = auxfLA;
    detail.fTI = auxfTI; detail.fD  = auxfD; detail.fEP = auxfEP;
    detail.fCV = auxfCV;
    fLI = rfLI; fLA = rfLA;
    fTI = rfTI; fD  = rfD; fEP = rfEP; fCV = rfCV;

end

function [mu, Sigma] = fixStats(X)
% X must be 2 x N

Xn = X';           % N x 2
mu = mean(Xn,1)';  % 2 x 1
Sigma = cov(Xn);   % 2 x 2

end
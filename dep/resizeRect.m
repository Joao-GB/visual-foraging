function newRect = resizeRect(rect, hFactor, vFactor)
% Mantém a posição do centro e multiplica os lados pelos fatores dados
    if nargin < 3, vFactor = hFactor; end
    isTransposed = false;
    if size(rect, 2) == 4
        isTransposed = true;
        rect = rect';
    end

    [cx, cy] = RectCenterd(rect);

    w = rect(3,:) - rect(1,:);
    h = rect(4,:) - rect(2,:);

    w = w * hFactor;
    h = h * vFactor;

    newRect = CenterRectOnPointd([zeros(1,numel(w)); zeros(1,numel(w)); w; h], cx, cy);
    if isTransposed, newRect = newRect'; end
end
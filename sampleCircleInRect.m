function center = sampleCircleInRect(ROIparams, totalRad)
    % ROIparams: [centerX, centerY, halfWidth, halfHeight]
    cx = ROIparams(1); cy = ROIparams(2);
    rx = ROIparams(3); ry = ROIparams(4);
    
    % 1. Calculate the "Safe Zone" dimensions
    safeWidth  = 2 * (rx - totalRad);
    safeHeight = 2 * (ry - totalRad);
    
    % 2. Check if the circle is physically too big for the rectangle
    if safeWidth < 0 || safeHeight < 0
        error('Erro: raio do círculo (%.2f) é maior que o retângulo permitido!', totalRad);
    end
    
    % 3. Sample uniformly in the safe zone
    % Start at (Left + Rad) and add a random portion of the SafeWidth
    x = (cx - rx + totalRad) + rand * safeWidth;
    y = (cy - ry + totalRad) + rand * safeHeight;
    
    center = [x, y];
end

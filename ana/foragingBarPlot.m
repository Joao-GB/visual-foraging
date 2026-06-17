function foragingBarPlot(s1,s2,plotIdx, x, y, N, col1, col2)
    L = numel(x);
    subplot(s1,s2,plotIdx);
    xc = categorical(string(x), string(x), 'Ordinal', true);
    b = bar(xc, [y;N-y]', 'stacked', 'FaceColor', 'flat');
    if ~isempty(col1), b(1).CData(col1(1),:) = col1(2:end); end
    b(2).CData = repmat(col2, L, 1);
end
L =length(fix.lims);
means = zeros(1,L);
for i = 1:L
  means(i) = mean(vel(fix.lims(1, i):fix.lims(2, i)));  
end
[what_1, where_1] = find(means > 16.5);
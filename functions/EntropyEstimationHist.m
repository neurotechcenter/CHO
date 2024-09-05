function E = EntropyEstimationHist(counts,centers)   
h = counts;
x = centers;

zero = find(h == 0 | isnan(h));
for k = 1:length(zero)
    h(zero(k)) = 1;
end
h = h/sum(h);
step = x(2)-x(1);
E = log(step)-sum(h.*log(h));
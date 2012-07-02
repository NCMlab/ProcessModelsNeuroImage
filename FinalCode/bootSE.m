function se = bootSE(data,N)

m = mean(data);
d = (data - m);
se = sqrt((d'*d)/length(data));
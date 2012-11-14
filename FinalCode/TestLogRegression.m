clear
N = 100;
M = randn(N,1);
X = sign(M);
a = 0.5;
M = M + a.*randn(N,1);
X = (X+1)./2;
 
% b = 0.01;
% c = 0.01;
% Y = c.*X + b.*M;
Nboot = 2000;
output = zeros(2,2,Nboot);
[b_logregPE se covb] = subfnLogisticRegress(X,M);
LOGt = b_logregPE./se
%b_regPE = subfnregress(X,M);
regSTATS = subfnregstats(X,M);
b_regPE = regSTATS.beta;
LRt = regSTATS.tstat.t


%corr([X M])
b_logreg = [0;0];
for i = 1:Nboot
    Samp =  floor(N*rand(N,1))+1;
    bootX = X(Samp);
    bootM = M(Samp);
    [b_logreg] = subfnLogisticRegress(bootX,bootM);
    b_reg = subfnregress(bootX,bootM);
    output(:,:,i) = [b_logreg b_reg];
end


figure(1)
clf
hist(squeeze(output(2,1,:)))
title(sprintf('Logistic regression, point est:%0.2f',b_logregPE(2)))

figure(2)
clf
hist(squeeze(output(2,2,:)))
title(sprintf('Linear regression, point est:%0.2f',b_regPE(2)))



function [S] = subfnLogisticRegressStats(Y,design)
% Perform logistic regression using the Newton-Raphson method
% A column of ones is added to the design.
% This function returns all statistics.
% 
% Written: 11/14/20012
% Written by: Jason Steffener
%
[N, Ncoef] = size(design);
X = [ones(N,1) design];
Ncoef = Ncoef + 1;
% Stopping criteria
tol = 10e-6;
% initial estimates 
sse = 10000 ;% using the sum of square errors for the stopping criteria
b0 = zeros(Ncoef,1);
p = ones(N,1).*0.5;
V = diag(0.5);
Niter = 0;
while sse > tol
    s = X'*(Y-p);
    m = (X'*V*X);
    beta = b0 + inv(m)*s;
    %beta = b0 + m\s; % This is supposed to be faste rthen using the inverse function
    p = exp(X*beta)./(1 + exp(X*beta));
    V = diag(p.*(1-p));
    err = (beta - b0);
    sse = err'*err;
    b0 = beta;
    Niter = Niter + 1;
end

covb = inv(X'*V*(X));
se = sqrt(diag(abs(covb)));
Z = beta./se;
dfe = N - Ncoef + 1
pval = 2*(tcdf(-abs(Z), 100000))
odds = exp(beta);
fit = X*beta;
% Classification
propCorClass = length(find(Y == round(p)))/N;
% Log likelihood
LL = sum(Y .* fit - log(1 + exp(fit)));
% Fit the reduced model
[beta0] = subfnLogisticRegress(Y,[]);
fit0 = ones(N,1)*beta0;
LL0 = sum(Y .* fit0 - log(1 + exp(fit0)));
LR = -2*(LL0 - LL);
% Likelihood ratio test
LRdf = length(beta) - 1;
LRp = 1 - chi2cdf(LR,LRdf);
McFadden = 1 - LL/LL0;
% Cox and Snell R^2
CoxSnell = 1 - (exp(LL0)/exp(LL))^(2/N);
% Nagelkerke R^2
Nagelkrk = (1 - (exp(LL0)/exp(LL))^(2/N) ) / (1 - exp(LL0)^(2/N));
S = {};
S.beta = beta;
S.covb = covb;
S.Zstat = {};
S.Zstat.Z = Z;
S.Zstat.se = se;
S.Zstat.odds = odds;
S.modelLL = LR;
S.modelLLdf = LRdf;
S.modelLLp = LRp;
S.m2LL = -2*LL;
S.McFadden = McFadden;
S.CoxSnell = CoxSnell;
S.Nagelkrk = Nagelkrk;
S.propCorClass = propCorClass;


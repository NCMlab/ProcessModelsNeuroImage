function beta = subfnLogisticRegress(Y,design)
% Perform logistic regression using the Newton-Raphson method
% A column of ones is added to the design.
% Only the beta values are estimated here.
% This function is written to just estimate the beta values for it is used
% with bootstrapping. For all stats use the subfnLogisticRegressStats
% function.
% 
% Written: 11/14/20012
% Written by: Jason Steffener
%
[~, Ncoef] = size(design);
N = length(Y);
X = [ones(N,1) design];
Ncoef = Ncoef + 1;
% Stopping criteria
tol = 10e-6;
% initial estimates 
sse = 10000 ;% using the sum of square errors for the stopping criteria
b0 = zeros(Ncoef,1);
p = ones(N,1).*0.5;
V = diag(0.5);
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
end


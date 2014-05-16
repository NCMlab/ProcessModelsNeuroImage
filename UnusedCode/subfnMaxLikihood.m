function [beta] = subfnMaxLikihood(Y,design)
% Perform logistic regression using the Newton-Raphson method
% Design needs a constant column
[N, Ncoef] = size(design);
X = [ones(N,1) design];
Ncoef = Ncoef + 1;
% Stopping criteria
tol = 10e-6;
% initial estimates 
sse = 10000 ;% using the sum of square errors for the stopping criteria
b0 = zeros(Ncoef,1);
p = zeros(N,1);

Niter = 0;
while sse > tol
    s = X'*(Y-p);
    m = (X'*X);
    beta = b0 + inv(m)*s;
    %beta = b0 + m\s; % This is supposed to be faste rthen using the inverse function
    p = X*beta;
    err = (beta - b0);
    sse = err'*err;
    b0 = beta;
    Niter = Niter + 1;
end
% 
% covb = inv(X'*V*(X));
% se = sqrt(diag(abs(covb)));
% fit = p;


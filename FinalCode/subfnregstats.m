function S = subfnregstats(y,design)  
% This is a version of the regstats MatLab program but slimmed down to only
% values needed.
%
% Add a column of ones to the design matrix
design = x2fx(design);
[Q,R] = qr(design,0);
  beta = R\(Q'*y);
  yhat = design*beta;
  residuals = y - yhat;
  nobs = length(y);
  p = length(beta);
  dfe = nobs-p;
  dft = nobs-1;
  %ybar = mean(y);
  sse = norm(residuals)^2;    % sum of squared errors
  %ssr = norm(yhat - ybar)^2;  % regression sum of squares
  %sst = norm(y - ybar)^2;     % total sum of squares;
  mse = sse./dfe;
  %h = sum(abs(Q).^2,2);
  %s_sqr_i = (dfe*mse - abs(residuals).^2./(1-h))./(dfe-1);
  %e_i = residuals./sqrt(s_sqr_i.*(1-h));
  ri = R\eye(p);
  xtxi = ri*ri';
  covb = xtxi*mse;
  se = sqrt(diag(covb));
  t = beta./se;
  pval = 2*(tcdf(-abs(t), dfe));
S = {};
S.beta = beta;
S.covb = covb;
S.tstat = {};
S.tstat.t = t;
S.tstat.se = se;
S.tstat.pval = pval;
S.tstat.dfe = dft;
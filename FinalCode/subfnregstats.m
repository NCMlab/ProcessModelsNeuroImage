function S = subfnregstats(y,design,Con,NoOnesFlag)
if nargin == 2
    Con = [];
    NoOnesFlag = 0;
elseif nargin ==3;
    NoOnesFlag = 0;
end
% This is a version of the regstats MatLab program but slimmed down to only
% values needed.
%
% Check to see if the y variable is dichotomous
if isequal(logical(y),y)
    S = subfnLogisticRegressStats(y,design);
else
    % Add a column of ones to the design matrix
    if ~NoOnesFlag 
        design = x2fx(design);
    end
    %design = [ones(size(design,1),1) design];
    [Q,R] = qr(design,0);
    beta = R\(Q'*y);
    yhat = design*beta;
    residuals = y - yhat;
    nobs = length(y);
    
    % create Hat matrix
    H = design*(R\Q');
    % calculate the leave one out cross validation statistic
    CV = ((residuals./(1-diag(H)))'*(residuals./(1-diag(H))))/nobs;
    p = length(beta);
    dfe = nobs-p;
    dft = nobs-1;
    ybar = mean(y);
    sse = norm(residuals)^2;    % sum of squared errors
    ssr = norm(yhat - ybar)^2;  % regression sum of squares
    sst = norm(y - ybar)^2;     % total sum of squares;
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
    S.B = beta.*(std(design)./std(y))';
    S.covb = covb;
    S.tstat = {};
    S.tstat.t = t;
    S.tstat.se = se;
    S.tstat.pval = pval;
    S.tstat.dfe = dft;
    S.rsquare = 1 - sse/sst;
    S.adjrsquare = 1 - (sse./sst)*(dft./dfe);
    S.fstat = {};
    S.fstat.sse = sse;
    S.fstat.dfe = dfe;
    S.fstat.dfr = p - 1;
    S.fstat.ssr = ssr;
    S.fstat.f = (ssr/S.fstat.dfr)/(sse/dfe);
    S.fstat.pval = fcdf(1/S.fstat.f, dfe, S.fstat.dfr);
    %S.AIC =  nobs*log(sum(residuals.^2)/nobs) + 2*p*(p+1)/(dfe-1) + 2*p;
    S.AIC = nobs*(log(2*pi) + 1 + log(sum(residuals.^2)/nobs)) + 2*(p + 1);
    S.BIC = nobs*(log(2*pi) + 1 + log(sum(residuals.^2)/nobs)) + (p + 1)*log(nobs);
    S.CV = CV;
    S.mse = mse;
    S.sse = sse;
    S.r = residuals;
    S.yhat = yhat;
end
if ~isempty(Con)
    [m n] = size(Con);
    if n ~= length(beta)
           error('Conrasts are not the same length as the beta vector.');
    else
        ConRes = {};
        for j = 1:m
            ConRes(j).c = Con(j,:);
            ConRes(j).con = (Con(j,:)*beta);
            ConRes(j).t = (Con(j,:)*beta)/sqrt(Con(j,:)*(covb)*Con(j,:)');
        end
    end
    S.Con = ConRes;
end

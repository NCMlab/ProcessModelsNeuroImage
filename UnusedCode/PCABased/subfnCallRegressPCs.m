function [beta, FVAL, EXITFLAG, fit] = subfnCallRegressPCs(data,coef)
if nargin == 1
    coef = [];
    coefFlag = 0;
else
    coefFlag = 1;
end
% This function makes the decision, based on the model, on how to estimate
% the parameters. The simple models can use linear regression but the
% models with interactions need to use an iterative process. 
% M: contains the SSFs of interest.
%
% define the options when using the fminsearch procedures 
options = optimset('TolX',1e-6);
%options = optimset('MaxIter',1000);
options = optimset(options,'MaxIter',1000);
options = optimset(options,'MaxFunEvals',1000);
options = optimset(options,'Algorithm','levenberg-marquardt');
options = optimset(options,'TolFun',1e-6);
options = optimset(options,'GradObj','on');
NSub = size(data.Y,1);
fit = {};
switch data.ModelNum
    case '4'
        % standard linear regression
        % Expected: 
        % M: [NSub x one or more SSFs]
        % X: [NSub x 1]
        % Y: [NSub x 1]
        % COV: [NSub x 0 or more]
        beta = subfnregress(data.Y,[data.M data.X data.COV]);
        FVAL = -9999;
        EXITFLAG = 1;
        fit = [ones(NSub,1) data.M data.X data.COV]*beta;
    case '14'
        NSSF = size(data.M,2);
        NMod = size(data.V,2);
        NCov = size(data.COV,2);
        if coefFlag == 0
            beta0 = zeros(1+NSSF+NMod+1+1+NCov,1);
        else
            beta0 = coef;
        end
        %beta0 = 1:1+NSSF+NMod+1+1+NCov
        [beta ,FVAL,EXITFLAG,OUTPUT] = fminsearch('subfnRegressPCs',beta0,'',data);
        
        NMed = NSSF;
        coef = beta;
        const = coef(1);        
        % coefficients for the mediator
        b = coef(1+1:1+NMed);
        % coefficient for the x effect
        cP = coef(1+NMed+1);        
        % coefficient(s) for the moderator
        v = coef(1+NMed+1+1:(NMod-1+NMed+1+2));
        % coefficients for the interaction term
        w = coef(1+NMed+1+NMod+1);
        % coefficients for the covariates
        b_cov = coef(end-NCov+1:end);
        if NCov > 0
            fit.Overall = ones(NSub,1)*const + data.M*b + data.X*cP + data.V*v + ((data.M*b).*(data.V*v))*w + data.COV*b_cov;
        else
            fit.Overall = ones(NSub,1)*const + data.M*b + data.X*cP + data.V*v + ((data.M*b).*(data.V*v))*w;
        end
        % Need to estimate the mediator and moderator and output them
        fit.mediator = data.M*b;
        fit.moderator = data.V*v;

    case '74'
        NSSF = size(data.M,2);
        NCov = size(data.COV,2);
        beta0 = zeros(1+NSSF+1+1+NCov,1);
        %beta0 = 1:1+NSSF+NMod+1+1+NCov
        [beta ,FVAL,EXITFLAG,OUTPUT] = fminsearch('subfnRegressPCs',beta0,'',data);
        
        NMed = NSSF;
        coef = beta;
        const = coef(1);        
        % coefficients for the mediator
        b = coef(1+1:1+NMed);
        % coefficient for the x effect
        cP = coef(1+NMed+1);        
        % coefficients for the interaction term
        w = coef(1+NMed+1+1);
        % coefficients for the covariates
        b_cov = coef(end-NCov+1:end);
        if NCov > 0
            fit.Overall = ones(N,1)*const + data.M*b + data.X*cP + ((data.M*b).*(data.X))*w + data.COV*b_cov;
        else
            fit.Overall = ones(N,1)*const + data.M*b + data.X*cP + ((data.M*b).*(data.X))*w;
        end
        % Need to estimate the mediator and moderator and output them
        fit.mediator = data.M*b;

end

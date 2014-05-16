function [fit, AIC] = subfnRegressPCs_ModelFit(coef,data)
ModelNum = data.ModelNum;
switch ModelNum
    case '4'
        % get data sizes
        N = length(data.Y);
        NMed = size(data.M,2);
        % coefficients for the mediator
        b = coef(1:NMed);
        % coefficient for the constant term in the model
        const = coef(NMed+1);
        % coefficient for the x effect
        cP = coef(NMed+2);
        % weighted sum of mediator PCs + X effect + weighted sum of moderator PCs +
        % interaction effect
        fit = data.M*b + data.X*cP + ones(N,1)*const;
        err = (data.Y - fit);
        F = err'*err;
        
    case '14'
        % get data sizes
        N = length(data.Y);
        NMed = size(data.M,2);
        NMod = size(data.V,2);
        NCov = size(data.COV,2);
        % coefficient for the constant term in the model
        const = coef(1);        
        % coefficients for the mediator
        b = coef(1+1:1+NMed);
        % coefficient for the x effect
        cP = coef(1+NMed+1);        
        % coefficient(s) for the moderator
        v = coef(1+NMed+1+1:(NMod-1+NMed+1+2))';
        % coefficients for the interaction term
        w = coef(1+NMed+1+NMod+1);
        % coefficients for the covariates
        b_cov = coef(end-NCov+1:end);

        % weighted sum of mediator PCs + X effect + weighted sum of moderator PCs +
        % interaction effect
        if NCov > 0
            fit = ones(N,1)*const + data.M*b + data.X*cP + data.V*v + ((data.M*b).*(data.V*v))*w + data.COV*b_cov;
        else
            fit = ones(N,1)*const + data.M*b' + data.X*cP + data.V*v' + ((data.M*b').*(data.V*v'))*w;
        end
        r = (data.Y - fit);
     
    p = length(coef);
    dfe = N - p;
    AIC =  N*log(sum(r.^2)/N) + 2*p*(p+1)/(dfe-1) + 2*p;
        
end
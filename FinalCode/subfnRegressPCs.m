function F = subfnRegressPCs(W,data)
ModelNum = data.ModelNum;
switch ModelNum
    case '4'
        % get data sizes
        N = length(data.Y);
        NMed = size(data.M,2);
        % coefficients for the mediator
        b = W(1:NMed);
        % coefficient for the constant term in the model
        const = W(NMed+1);
        % coefficient for the x effect
        cP = W(NMed+2);
        
        Y = data.Y;
        M = data.M;
        X = data.X;
        % weighted sum of mediator PCs + X effect + weighted sum of moderator PCs +
        % interaction effect
        fit = M*b + X*cP + ones(N,1)*const;
        err = (Y - fit);
        F = err'*err;
    case '7'
        % get data sizes
        N = length(data.Y);
        NMed = size(data.M,2);
        NMod = size(data.V,2);
        % coefficients for the mediator
        b = W(1:NMed);
        % coefficient(s) for the moderator
        v = W(NMed+1:NMed+NMod);
        % coefficients for the interaction term
        w = W(NMed+NMod+1);
        % coefficient for the constant term in the model
        const = W(NMed+NMod+2);
        % coefficient for the x effect
        cP = W(NMed+NMod+3);
        
        Y = data.Y;
        M = data.M;
        V = data.V;
        X = data.X;
        % weighted sum of mediator PCs + X effect + weighted sum of moderator PCs +
        % interaction effect
        fit = M*b + X*cP + V*v + ((M*b).*(V*v))*w + ones(N,1)*const;
        err = (Y - fit);
        F = err'*err;
end
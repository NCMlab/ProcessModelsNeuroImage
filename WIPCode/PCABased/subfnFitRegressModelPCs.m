function fit = subfnFitRegressModelPCs(data,beta,ModelNum)
% This function makes the decision, based on the model, on how to fit the
% model using precalculated beta weights
% M: contains the SSFs of interest.
%


NSub = size(data.Y,1);
switch ModelNum
    case '4'
        % standard linear regression
        % Expected: 
        % M: [NSub x one or more SSFs]
        % X: [NSub x 1]
        % Y: [NSub x 1]
        % COV: [NSub x 0 or more]

        fit = [ones(NSub,1) data.M data.X data.COV]*beta;
    case '14'
        NSSF = size(data.M,2);
        NMod = size(data.V,2);
        NCov = size(data.COV,2);
        beta0 = zeros(1+NSSF+NMod+1+1+NCov,1);
        beta0 = 1:1+NSSF+NMod+1+1+NCov
%        [beta ,FVAL,EXITFLAG,OUTPUT] = fminsearch('subfnRegressPCs',beta0,options,data);
        %fit = ones(Nsub,1)*const + data.M*b + data.X*cP + data.V*v + ((data.M*b).*(data.V*v))*w + data.COV;
end

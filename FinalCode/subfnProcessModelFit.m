function [ParameterToBS Parameters] = subfnProcessModelFit(data,ModelNum,PointEst)
ParameterToBS =[];
Parameters = {};
switch ModelNum
    case '4'
        % This is the simple mediation case which can handle covariates on
        % M and Y and multiple mediators, M.
        Nmed = size(data.M,2);
        Ndata = size(data.Y,1);
        a = zeros(Nmed,1);
        % covariates
        if size(data.COV,2) > 0 
            for i = 1:Nmed
                temp1 = regress(data.M(:,i),[data.X data.COV ones(Ndata,1)]);
                a(i) = temp1(1);
            end
            temp2 = regress(data.Y,[data.M data.X data.COV ones(Ndata,1)]);
            b = temp2(1:Nmed);
        % no covariates    
        else 
            for i = 1:Nmed
                temp1 = regress(data.M(:,i),[data.X ones(Ndata,1)]);
                a(i) = temp1(1);
            end
            temp2 = regress(data.Y,[data.M data.X ones(Ndata,1)]);
            b = temp2(1:Nmed);
        end
        % the indirect effect which will be bootstrapped
        ab = a.*b;
        ParameterToBS = ab;
        % Now all parameters of interest for the model are calculated.
        if PointEst 
            S1 = cell(Nmed,1);
            if size(data.COV,2) > 0 % covariates
                S2 = subfnregstats(data.Y,[data.M data.X data.COV]);
                for i = 1:Nmed
                    S1{i} = subfnregstats(data.M(:,i),[data.X data.COV]);
                end
                S3 = subfnregstats(data.Y,[data.X data.COV]);
            else
                S2 = subfnregstats(data.Y,[data.M data.X]);
                 for i = 1:Nmed
                    S1{i} = subfnregstats(data.M(:,i),[data.X ]);
                end
                S3 = subfnregstats(data.Y,data.X);    
            end
            Parameters = {};
            Parameters.A = {};
            Parameters.B = {};
            Parameters.AB = {};
            Parameters.C = {};
            Parameters.CP = {};
            for i = 1:Nmed
                Parameters.A{i}.beta = S1{i}.beta(2);
                Parameters.A{i}.se = S1{i}.tstat.se(2);
                Parameters.A{i}.t = S1{i}.tstat.t(2);
                Parameters.A{i}.p = S1{i}.tstat.pval(2);
                Parameters.A{i}.df = S1{i}.tstat.dfe;
                Parameters.B{i}.beta = S2.beta(i + 1);
                Parameters.B{i}.se = S2.tstat.se(i + 1);
                Parameters.B{i}.t = S2.tstat.t(i + 1);
                Parameters.B{i}.p = S2.tstat.pval(i + 1);
                Parameters.B{i}.df = S2.tstat.dfe;
                Parameters.AB{i}.pointEst = ab(i);
            end
            Parameters.C.beta = S3.beta(2);
            Parameters.C.se = S3.tstat.se(2);
            Parameters.C.t = S3.tstat.t(2);
            Parameters.C.p = S3.tstat.pval(2);
            Parameters.C.df = S3.tstat.dfe;
            Parameters.CP.beta = S2.beta(end);
            Parameters.CP.se = S2.tstat.se(end);
            Parameters.CP.t = S2.tstat.t(end);
            Parameters.CP.p = S2.tstat.pval(end);
            Parameters.CP.df = S2.tstat.dfe;
        end
        
    case '14'
        % This is the moderated mediation model. The moderation (V) occurs
        % between the mediator(s) M and Y. The conditional effect of X on Y
        % via M is evaluated at multiple moderation values. The confidence
        % intervals for each of these moderating values are calculated via
        % bootstrapping. Since the model is calulated 22 (21 steps + 1 point 
        % estimate) times this model is then 22 times slower than the
        % simple mediation model.
        Nsteps = 11;
        Nmed = size(data.M,2);
        Ndata = size(data.Y,1);
        a = zeros(Nmed,1);
        % Find the values of the moderator for probing
        % what needs to be bootstrapped here is the point estimate but also
        % 20 values of the moderator and the J-N values.
        minV = min(data.V);
        maxV = max(data.V);
        rangeV = maxV - minV;
        stepV = rangeV/(Nsteps -1);
        probeV = [0 minV:stepV:maxV];
        ParameterToBS = zeros(Nmed,Nsteps);
        if size(data.COV,2) > 0 % covariates
            for k = 1:Nsteps + 1
                Interaction = zeros(Ndata,Nmed);
                for j = 1:Nmed
                    Interaction(:,j) = data.M(:,j).*(data.V - probeV(k));
                    S1 = subfnregstats(data.M(:,j),[data.X data.COV]);
                    a(j) = S1.beta(2);
                end
                S2 = subfnregstats(data.Y,[data.M (data.V - probeV(k)) Interaction data.X data.COV]);
                ParameterToBS(:,k) = a.*(S2.beta(2:Nmed+1));
            end
        else % no covariates
            for k = 1:Nsteps + 1
                Interaction = zeros(Ndata,Nmed);
                for j = 1:Nmed
                    Interaction(:,j) = data.M(:,j).*(data.V - probeV(k));
                    S1 = subfnregstats(data.M(:,j),[data.X]);
                    a(j) = S1.beta(2);
                end
                S2 = subfnregstats(data.Y,[data.M (data.V - probeV(k)) Interaction data.X]);
                ParameterToBS(:,k) = a.*(S2.beta(2:Nmed+1));
            end
        end
        Parameters = {};
        if PointEst
            S1 = {};
            if size(data.COV,2) > 0 % covariates
                Interaction = zeros(Ndata,Nmed);
                for j = 1:Nmed
                    Interaction(:,j) = data.M(:,j).*(data.V);
                    S1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
                end
                S2 = subfnregstats(data.Y,[data.M (data.V) Interaction data.X data.COV]);
                S3 = subfnregstats(data.Y,[data.X data.COV]);
            else % no covariates
                Interaction = zeros(Ndata,Nmed);
                for j = 1:Nmed
                    Interaction(:,j) = data.M(:,j).*(data.V);
                    S1{j} = subfnregstats(data.M(:,j),[data.X]);
                end
                S2 = subfnregstats(data.Y,[data.M (data.V) Interaction data.X]);
                S3 = subfnregstats(data.Y,[data.X]);
            end
            
            Parameters.A = {};
            Parameters.B = {};
            Parameters.V = {};
            Parameters.Int = {};
            Parameters.AB = {};
            Parameters.C = {};
            Parameters.CP = {};
            for i = 1:Nmed
                Parameters.A{i}.beta = S1{i}.beta(2);
                Parameters.A{i}.se = S1{i}.tstat.se(2);
                Parameters.A{i}.t = S1{i}.tstat.t(2);
                Parameters.A{i}.p = S1{i}.tstat.pval(2);
                Parameters.A{i}.df = S1{i}.tstat.dfe;
                Parameters.B{i}.beta = S2.beta(i + 1);
                Parameters.B{i}.se = S2.tstat.se(i + 1);
                Parameters.B{i}.t = S2.tstat.t(i + 1);
                Parameters.B{i}.p = S2.tstat.pval(i + 1);
                Parameters.B{i}.df = S2.tstat.dfe;
                
                Parameters.Int{i}.beta = S2.beta(1 + Nmed + 1 + i);
                Parameters.Int{i}.se = S2.tstat.se(1 + Nmed + 1 + i);
                Parameters.Int{i}.t = S2.tstat.t(1 + Nmed + 1 + i);
                Parameters.Int{i}.p = S2.tstat.pval(1 + Nmed + 1 + i);
                Parameters.Int{i}.df = S2.tstat.dfe;
                
                for k = 1:length(probeV)
                    Parameters.AB{k,i}.pointEst = ParameterToBS(i,k);
                    Parameters.AB{k,i}.V = probeV(k);
                end
            end
            Parameters.V{1}.beta = S2.beta(Nmed+2);
            Parameters.V{1}.se = S2.tstat.se(Nmed+2);
            Parameters.V{1}.t = S2.tstat.t(Nmed+2);
            Parameters.V{1}.p = S2.tstat.pval(Nmed+2);
            Parameters.V{1}.df = S2.tstat.dfe;

            Parameters.C.beta = S3.beta(2);
            Parameters.C.se = S3.tstat.se(2);
            Parameters.C.t = S3.tstat.t(2);
            Parameters.C.p = S3.tstat.pval(2);
            Parameters.C.df = S3.tstat.dfe;

            Parameters.CP.beta = S2.beta(1+Nmed+1+Nmed+1);
            Parameters.CP.se = S2.tstat.se(1+Nmed+1+Nmed+1);
            Parameters.CP.t = S2.tstat.t(1+Nmed+1+Nmed+1);
            Parameters.CP.p = S2.tstat.pval(1+Nmed+1+Nmed+1);
            Parameters.CP.df = S2.tstat.dfe;
        end
            
end
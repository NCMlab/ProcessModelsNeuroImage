function [ParameterToBS Parameters] = subfnProcessModelFit(data,PointEstFlag)
% PointEstFlag:
%           set to 0 if you only want to calculate the point estimate for
%           the model. This is used for each bootstrapping test.
%           set to 1 to estimate all parameters in the model.
%


Parameters = {};
Nsteps = 11;
switch data.ModelNum
    case '1'
        ParameterToBS = struct('names','CondMod','values',zeros(1,Nsteps + 1),'probeValues',zeros(1,Nsteps + 1),'probeMod',0);
       
        Ndata = size(data.Y,1);
        % the code below works if there is a covariate or not
        
        % whether or not to run the regression at multiple values of the moderator
        % First, check to see if the interaction effect is significant or
        % not.
        Model1 = subfnregstats(data.Y,[data.X data.M (data.M).*data.X data.COV]);
        % When this program is called during boot strapping it needs to
        % know whether or not to probe the interaction.
        % It should only check to see if the interaction is significant for
        % when the point estimate is being tested and not for any boot
        % strap re-estimates.
        if Model1.tstat.pval(4) < max(data.Thresholds)
            ParameterToBS.probeMod = 1;
        end
        %
        % If the probeMd flag is set to TRUE then check to see if the
        % interaction is significant. If so then leave the probeMod flag
        % set to TRUE. If the interaction is not significant then set the
        % flag to FALSE.
        if data.ProbeMod %
            % check to see if the interaction is significant! This should
            % only be checked the first time through.
            % if S.tstat.pval(4) < max(data.Thresholds)
            ParameterToBS.values(1,1) = Model1.beta(2);
            minM = min(data.M);
            maxM = max(data.M);
            rangeM = maxM - minM;
            stepM = rangeM/(Nsteps -1);
            probeM = [0 minM:stepM:maxM];
            for j = 2:Nsteps + 1
                temp = subfnregress(data.Y,[data.X (data.M-probeM(j))  (data.M-probeM(j)).*data.X data.COV]);
                ParameterToBS.values(1,j) = temp(2);
                ParameterToBS.probeValues(1,j) = probeM(j);
            end
        else
            ParameterToBS.values = Model1.beta(2);
            ParameterToBS.probeValues = 0;
        end
        
        % Now all parameters of interest for the model are calculated.
        if PointEstFlag
            % also fit the direct model without the modulator in it
            % Model with NO Interaction
            Model2 = subfnregstats(data.Y,[data.X data.M data.COV]);
            % Find the R2 increase due to the interaction
            diffS = subfnCalculateModelFitDiff(Model1,Model2);
            Parameters = {};
            Parameters.Model1.const = subfnSetParameters('const', Model1, 1);
            Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,2);',data.Xname,data.Xname);
            eval(Str)
            Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,3);',data.Mname,data.Mname);
            eval(Str)
            Str = sprintf('Parameters.Model1.%s_x_%s=subfnSetParameters(''%s_x%s'',Model1,4);',data.Xname,data.Mname,data.Xname,data.Mname);
            eval(Str)
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,4+j);',data.COVname{j},data.COVname{j});
                eval(Str)
            end
            Parameters.Model1.Outcome = data.Yname;
            Parameters.Model1.Model = subfnSetModelParameters(Model1);
            % Model 2, no interaction
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,2);',data.Xname,data.Xname);
            eval(Str)
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,3);',data.Mname,data.Mname);
            eval(Str)
            
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,3+j);',data.COVname{j},data.COVname{j});
                eval(Str)
            end
            Parameters.Model2.Outcome = data.Yname;
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            Parameters.DiffModel = subfnSetModelParameters(diffS);
            tcrit = tinv(1 - max(data.Thresholds)/2,length(data.X) - (4 + size(data.COV,2)));
            Parameters.JNvalue = subfnJohnsonNeyman(Model1.beta(2),Model1.covb(2,2),Model1.beta(4),Model1.covb(4,4),Model1.covb(2,4),tcrit);
        end

    case '4'
        % This is the simple mediation case which can handle covariates on
        % M and Y and multiple mediators, M.
        
        Nmed = size(data.M,2);
        Ndata = size(data.Y,1);
        NameStruct = cell(Nmed,1);
        for j = 1:Nmed 
            NameStruct{j} = sprintf('AB%d',j);
        end
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,Nsteps + 1),'probeValues',zeros(1,1),'probeMod',0);
        % redo the setup so that results are organized based on models
        % =================================================================
        % Model1
        a = zeros(Nmed,1);
        for i = 1:Nmed
            % A branch model
            tempModel1 = subfnregress(data.M(:,i),[data.X data.COV]);
            a(i) = tempModel1(2);
        end
        % =================================================================
        % Model2
        % B branch model
        tempModel2 = subfnregress(data.Y,[data.M data.X data.COV]);
        b = tempModel2(2:Nmed+1);
        % =================================================================
        % Bootstrap values
        % the indirect effect which will be bootstrapped
        ab = a.*b;
        ParameterToBS.values = ab;
        ParameterToBS.k2 = subfnCalculateKappa2(data.X, data.M, data.Y, a, b);
        % Now all parameters of interest for the model are calculated.
        if PointEstFlag 
            S1 = cell(Nmed,1);
            % B branch model
            Model2 = subfnregstats(data.Y,[data.M data.X data.COV]);
            for i = 1:Nmed
                % A branch model
                Model1{i} = subfnregstats(data.M(:,i),[data.X data.COV]);
            end
            % C branch model
            Model3 = subfnregstats(data.Y,[data.X data.COV]);
            % Fill in the Parameters structure with all results 
            % from Model 1
            Parameters = {};
            for i = 1:Nmed
                Parameters.Model1{i}.const = subfnSetParameters('const', Model1{i}, 1);
                Str = sprintf('Parameters.Model1{i}.%s=subfnSetParameters(''%s'',Model1{i},2);',data.Xname,data.Xname);
                eval(Str)
                for j = 1:size(data.COV,2)
                    Str = sprintf('Parameters.Model1{i}.%s=subfnSetParameters(''%s'',Model1{i},2+j);',data.COVname{j},data.COVname{j});
                    eval(Str)
                end
                Parameters.Model1{i}.Model = subfnSetModelParameters(Model1{i});
            end
            Parameters.Model1{i}.Outcome = data.Mname;

            % Fill in the Parameters structure with all results
            % from Model 2
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            for j = 1:Nmed
                Str = sprintf('Parameters.Model2.%s%d=subfnSetParameters(''%s'',Model2,1+j);',data.Mname,j,data.Mname);
                eval(Str);
            end
            Str = sprintf('Parameters.Model2.%s = subfnSetParameters(''%s'', Model2, 1+Nmed+1);',data.Xname,data.Xname);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+1+j);',data.COVname{j},data.COVname{j});
                eval(Str);
            end
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            Parameters.Model2.Outcome = data.Yname;
            % Fill in the Parameters structure with all results
            % from Model 3            
            Parameters.Model3.const = subfnSetParameters('const', Model3, 1);
            Str = sprintf('Parameters.Model3.%s = subfnSetParameters(''%s'', Model3, 2);',data.Xname,data.Xname);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+1+j);',data.COVname{j},data.COVname{j});
                eval(Str);
            end
            Parameters.Model3.Model = subfnSetModelParameters(Model3);
            % Fill in the Parameters structure with all results
            % from the bootstrapped values
            for i = 1:Nmed
                Str = sprintf('Parameters.%s.pointEst = ParameterToBS.values(i);',ParameterToBS.names(i,:));
                eval(Str);
            end
            Parameters.JohnsonNeyman = -99;
            Parameters.Model3.Outcome = data.Yname;

        end
        
    case '7'
        % this is the model where the branch from X to M is moderated
        [Ndata Nmed] = size(data.M);
        NameStruct = cell(Nmed,1);
        for j = 1:Nmed 
            NameStruct{j} = sprintf('CondAB%d',j);
        end
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,Nsteps + 1),'probeValues',zeros(1,Nsteps + 1),'probeMod',0);

        a = zeros(Nmed,1);
        % create the interaction terms
        S1 = cell(Nmed,1);
        for j = 1:Nmed
            Interaction = data.M(:,j).*(data.W);
            S1{j} = subfnregstats(data.M(:,j),[data.X Interaction data.W data.COV]);
            % check each interaction term and set the flag here for whether
            % to probe the interaction.
            if S1{j}.tstat.pval(3) < max(data.Thresholds)
                ParameterToBS.probeMod = 1;
            end
        end
        % B branch model
        S2 = subfnregress(data.Y,[data.M data.X data.COV]);
        if data.ProbeMod %
            % check to see if the interaction is significant! This should
            % only be checked the first time through.
            % if S.tstat.pval(4) < max(data.Thresholds)
            for j = 1:Nmed
                ParameterToBS.values(j,1) = S1{j}.beta(2)*S2(2);
            end
            minW = min(data.W);
            maxW = max(data.W);
            rangeW = maxW - minW;
            stepW = rangeW/(Nsteps -1);
            probeW = [0 minW:stepW:maxW];
            for k = 2:Nsteps + 1
                for j = 1:Nmed
                    Interaction = data.M(:,j).*(data.W - probeW(k));
                    temp = subfnregress(data.M(:,j),[data.X Interaction (data.W-probeW(k)) data.COV]);
                    ParameterToBS.values(j,k) = temp(2)*S2(2);
                    ParameterToBS.probeValues(1,k) = probeW(k);
                end
            end
        else
            temp = zeros(Nmed,1);
            for j = 1:Nmed
                temp(j) = S1{j}.beta(2)*S2(2);
            end
            ParameterToBS.values = temp;
            ParameterToBS.probeValues = 0;
        end
        Parameters = {};
        if PointEstFlag
            S2 = subfnregstats(data.Y,[data.M data.X data.COV]);
            S3 = subfnregstats(data.Y,[data.X data.COV]);
            %noIntS = subfnregstats(data.Y,[data.M data.V  data.X data.COV]);
            %diffS = subfnCalculateModelFitDiff(S,noIntS);
            %Parameters.EffOfInt = subfnSetModelParameters(diffS);
            Parameters.CModel = subfnSetModelParameters(S3);
            Parameters.BModel = subfnSetModelParameters(S2);
            % calculate the Johnson-Neyman value
%            JNvalue = subfnJohnsonNeyman(S.beta(2),S.covb(2,2),S.beta(4),S.covb(4,4),S.covb(2,4),data.tcrit);
%             Parameters.JohnsonNeyman = JNvalue;
            for j = 1:Nmed
                str = sprintf('Parameters.A%dModel = subfnSetModelParameters(S1{j});',j); eval(str);
               
                str = sprintf('Parameters.A%d = subfnSetParameters(''A%d'',S1{j},2);',j,j);eval(str);
                str = sprintf('Parameters.A%dconst = subfnSetParameters(''A%dconst'',S1{j},1);',j,j);eval(str);
                str = sprintf('Parameters.Int%d = subfnSetParameters(''Int%d'',S1{j},3);',j,j);eval(str);
                str = sprintf('Parameters.W%d = subfnSetParameters(''W%d'',S1{j},4);',j,j);eval(str);
                str = sprintf('Parameters.B%d = subfnSetParameters(''V%d'',S2,1+j);',j,j);eval(str);
            end
             Parameters.Bconst = subfnSetParameters('Bconst',S2,1);
             Parameters.CP = subfnSetParameters('CP',S2,Nmed+2);
             Parameters.C = subfnSetParameters('C',S3,2);
             Parameters.Cconst = subfnSetParameters('Cconst',S3,1);
        end
    case '14'
        % This is the moderated mediation model. The moderation (V) occurs
        % between the mediator(s) M and Y. The conditional effect of X on Y
        % via M is evaluated at multiple moderation values. The confidence
        % intervals for each of these moderating values are calculated via
        % bootstrapping. Since the model is calulated 22 (21 steps + 1 point 
        % estimate) times this model is then 22 times slower than the
        % simple mediation model.

        [Ndata Nmed] = size(data.M);
        NameStruct = cell(Nmed,1);
        for j = 1:Nmed 
            NameStruct{j} = sprintf('CondAB%d',j);
        end
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,Nsteps + 1),'probeValues',zeros(1,Nsteps + 1),'probeMod',0);
        a = zeros(Nmed,1);
        % First, check to see if the interaction effect is significant or
        % not.
        Interaction = zeros(Ndata,Nmed);
        S1 = {};
        for j = 1:Nmed
            Interaction(:,j) = data.M(:,j).*(data.V);
            S1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
            a(j) = S1{j}.beta(2);
        end
        S = subfnregstats(data.Y,[data.M data.V  Interaction data.X data.COV]);
        % check to see if the interaction is significant
        
        for j = 1:Nmed
            if S.tstat.pval(1+Nmed+1+j) < max(data.Thresholds)
                ParameterToBS.probeMod = 1;
            end
        end
        if data.ProbeMod
            ParameterToBS.values(1,1) = a.*(S.beta(2:Nmed+1));
            minV = min(data.V);
            maxV = max(data.V);
            rangeV = maxV - minV;
            stepV = rangeV/(Nsteps -1);
            probeV = [0 minV:stepV:maxV];
            for k = 2:Nsteps + 1
                Interaction = zeros(Ndata,Nmed);
                for j = 1:Nmed
                    Interaction(:,j) = data.M(:,j).*(data.V - probeV(k));
                    S1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
                    a(j) = S1{j}.beta(2);
                end
                 % B branch model
                S2 = subfnregstats(data.Y,[data.M (data.V - probeV(k)) Interaction data.X data.COV]);
                ParameterToBS.values(:,k) = a.*(S2.beta(2:Nmed+1));
            end
            ParameterToBS.probeValues = probeV;
        else
            ParameterToBS.values = a.*(S.beta(2:Nmed+1));
            ParameterToBS.probeValues = 0;
        end
        Parameters = {};
        if PointEstFlag
            S3 = subfnregstats(data.Y,[data.X data.COV]);
            noIntS = subfnregstats(data.Y,[data.M data.V  data.X data.COV]);
            diffS = subfnCalculateModelFitDiff(S,noIntS);
            Parameters.EffOfInt = subfnSetModelParameters(diffS);
            Parameters.CModel = subfnSetModelParameters(S3);
            Parameters.BModel = subfnSetModelParameters(S);
            % calculate the Johnson-Neyman value
            JNvalue = subfnJohnsonNeyman(S.beta(2),S.covb(2,2),S.beta(4),S.covb(4,4),S.covb(2,4),data.tcrit);
            Parameters.JohnsonNeyman = JNvalue;
            for i = 1:Nmed
                Parameters.AModel = subfnSetModelParameters(S1{i});
                Parameters.A{i} = subfnSetParameters('A', S1{i},2);
                Parameters.Aconst{i} = subfnSetParameters('Aconst', S1{i},1);
                Parameters.B{i} = subfnSetParameters('B', S,i + 1);
                Parameters.Int{i} = subfnSetParameters('Int',S,1 + Nmed + 1 + i);
            end
            Parameters.Bconst = subfnSetParameters('Bconst',S,1);
            Parameters.V = subfnSetParameters('V',S,Nmed+2);
            Parameters.C = subfnSetParameters('C',S3,2);
            Parameters.Cconst = subfnSetParameters('Cconst',S3,1);
            Parameters.CP = subfnSetParameters('Cconst',S,1+Nmed+1+Nmed+1);
            
        end
end

if PointEstFlag
    Parameters.Xname = data.Xname;
    Parameters.Mname = data.Mname;
    Parameters.Yname = data.Yname;
    Parameters.Vname = data.Vname;
    Parameters.Wname = data.Wname;
    Parameters.ModelNum = data.ModelNum;
    Parameters.SampleSize = length(data.X);
end
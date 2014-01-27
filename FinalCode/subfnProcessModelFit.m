function [ParameterToBS Parameters] = subfnProcessModelFit(data,PointEstFlag)
% PointEstFlag:
%           set to 0 if you only want to calculate the point estimate for
%           the model. This is used for each bootstrapping test.
%           set to 1 to estimate all parameters in the model.
%
% MODEL 1 is the A Branch. This is from X to M
% MODEL 2 is the B Branch. This is from M to Y
% Model 3 is the C Branch. This is from X to Y

Parameters = {};

switch data.ModelNum
    case '1'
        %
        %     M
        %     |
        % X ----- Y
        %
        % Check to see if the moderator, in this case M, is dicotomous.
        minM = min(data.M);
        maxM = max(data.M);
        rangeM = maxM - minM;
        if rangeM == 1
            probeM = [0 1];
        else           
            probeM = subfnCalculateProbeValues(data.M);
        end

        ParameterToBS = struct('names','CondMod','values',zeros(1,length(probeM)),'probeValues',zeros(1,length(probeM)),'ProbeMod',0);
        Ndata = size(data.Y,1);
        % First, check to see if the interaction effect is significant or
        % not.
        Model1 = subfnregstats(data.Y,[data.X data.M (data.M).*data.X data.COV]);
        % When this program is called it needs to
        % know whether or not to probe the interaction.
        % It should only check to see if the interaction is significant for
        % when the point estimate is being tested and not for any boot
        % strap re-estimates.

        ParameterToBS.ProbeMod = 1;
        %
        % If the probeMd flag is set to TRUE then check to see if the
        % interaction is significant. If so then leave the ProbeMod flag
        % set to TRUE. If the interaction is not significant then set the
        % flag to FALSE.
        if data.ProbeMod %
            % check to see if the interaction is significant! This should
            % only be checked the first time through.
            % if S.tstat.pval(4) < max(data.Thresholds)
            
            ParameterToBS.values(1,1) = Model1.beta(2);
            
            
            
            for j = 1:length(probeM)
                %temp = subfnregress(data.Y,[data.X (data.M-probeM(j))  (data.M-probeM(j)).*data.X data.COV]);
                ParameterToBS.values(1,j) = Model1.beta(2)+Model1.beta(4)*probeM(j);
                ParameterToBS.probeValues(1,j) = probeM(j);
            end
        else
            ParameterToBS.values = Model1.beta(2);
            ParameterToBS.probeValues = 0;
        end
        ParameterToBS.k2 = 0;
        % Now all parameters of interest for the model are calculated.
        if PointEstFlag
            % also fit the direct model without the modulator in it

            % Model with NO Interaction
            Model2 = subfnregstats(data.Y,[data.X data.M data.COV]);
            % Find the R2 increase due to the interaction
            diffS = subfnCalculateModelFitDiff(Model1,Model2);
            Parameters = {};
            Parameters.Model1.const = subfnSetParameters('const', Model1, 1);
            Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,2);',data.names.X,data.names.X);
            eval(Str)
            Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,3);',data.names.M{1},data.names.M{1});
            eval(Str)
            Str = sprintf('Parameters.Model1.%s_x_%s=subfnSetParameters(''%s_x%s'',Model1,4);',data.names.X,data.names.M{1},data.names.X,data.names.M{1});
            eval(Str)
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,4+j);',data.names.COV{j},data.names.COV{j});
                eval(Str)
            end
            Parameters.Model1.Outcome = data.names.Y;
            Parameters.Model1.Model = subfnSetModelParameters(Model1);
            % Model 2, no interaction
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,2);',data.names.X,data.names.X);
            eval(Str)
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,3);',data.names.M{1},data.names.M{1});
            eval(Str)
            
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,3+j);',data.names.COV{j},data.names.COV{j});
                eval(Str)
            end
            Parameters.Model2.Outcome = data.names.Y;
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            Parameters.DiffModel = subfnSetModelParameters(diffS);
            tcrit = tinv(1 - max(data.Thresholds)/2,length(data.X) - (4 + size(data.COV,2)));
            Parameters.JNvalue = subfnJohnsonNeyman(Model1.beta(2),Model1.covb(2,2),Model1.beta(4),Model1.covb(4,4),Model1.covb(2,4),tcrit);

        end
        
    case '4'
        %
        %     M
        %    / \
        %   /   \
        %  X     Y
        %
        % This is the simple mediation case which can handle covariates on
        % M and Y and multiple mediators, M.
        
        Nmed = size(data.M,2);
        Ndata = size(data.Y,1);
        NameStruct = cell(Nmed,1);
        for j = 1:Nmed
            NameStruct{j} = sprintf('AB%d',j);
        end
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,1),'probeValues',zeros(1,1),'ProbeMod',0);
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
                Str = sprintf('Parameters.Model1{i}.%s=subfnSetParameters(''%s'',Model1{i},2);',data.names.X,data.names.X);
                eval(Str)
                for j = 1:size(data.COV,2)
                    Str = sprintf('Parameters.Model1{i}.%s=subfnSetParameters(''%s'',Model1{i},2+j);',data.names.COV{j},data.names.COV{j});
                    eval(Str)
                end
                Parameters.Model1{i}.Model = subfnSetModelParameters(Model1{i});
                Parameters.Model1{i}.Outcome = data.names.M{1};
            end
            
            
            % Fill in the Parameters structure with all results
            % from Model 2
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            for j = 1:Nmed
                Str = sprintf('Parameters.Model2.%s%d=subfnSetParameters(''%s'',Model2,1+j);',data.names.M{j},j,data.names.M{j});
                eval(Str);
            end
            Str = sprintf('Parameters.Model2.%s = subfnSetParameters(''%s'', Model2, 1+Nmed+1);',data.names.X,data.names.X);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+1+j);',data.names.COV{j},data.names.COV{j});
                eval(Str);
            end
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            Parameters.Model2.Outcome = data.names.Y;
            % Fill in the Parameters structure with all results
            % from Model 3
            Parameters.Model3.const = subfnSetParameters('const', Model3, 1);
            Str = sprintf('Parameters.Model3.%s = subfnSetParameters(''%s'', Model3, 2);',data.names.X,data.names.X);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+1+j);',data.names.COV{j},data.names.COV{j});
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
            Parameters.Model3.Outcome = data.names.Y;
            
        end
    case '6'
       %    M1--M2
       %   /      \
       %  /        \
       % X          Y
       %
%        Model 1: X to M1
%        Model 2: X and M1 to M2
%        Model 3: X, M1 and M2 to Y
%        Model 4: X to Y
%        Only the full indirect path is being testing with bootstraping. I
%        am not sure if I should put in the shorter paths as well. They are
%        actually calculable with tests of Model 4.

        Nmed = size(data.M,2);
        Ndata = size(data.Y,1);

        NameStruct = {'M1' 'M1M2'};

        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,1),'probeValues',zeros(1,1),'ProbeMod',0);
        % redo the setup so that results are organized based on models
        % =================================================================
        % Model 1
        tempModel1 = subfnregress(data.M(:,1),[data.X data.COV]);
        % =================================================================
        % Model 2
        tempModel2 = subfnregress(data.M(:,2),[data.M(:,1) data.X data.COV]);
        % =================================================================
        % Model 3
        tempModel3 = subfnregress(data.Y,[data.M(:,2) data.M(:,1) data.X data.COV]);
        % =================================================================
        % Model 4
        tempModel4 = subfnregress(data.Y,[data.X data.COV]);

        % =================================================================
        % Bootstrap values
        % the indirect effect which will be bootstrapped
        % X - M1 - M2 - Y
        ParameterToBS.values(1) = tempModel1(2)*tempModel2(2)*tempModel3(2);
        ParameterToBS.values(2) = tempModel1(2)*tempModel3(2);
        ParameterToBS.k2 = 0;
        if PointEstFlag
            Model1 = subfnregstats(data.M(:,1),[data.X data.COV]);
            Model2 = subfnregstats(data.M(:,2),[data.M(:,1) data.X data.COV]);
            Model3 = subfnregstats(data.Y,[data.M(:,2) data.M(:,1) data.X data.COV]);
            Model4 = subfnregstats(data.Y,[data.X data.COV]);
            Parameters = {};

            % Fill in the Parameters structure with all results
            % Model 1
            Parameters.Model1.const = subfnSetParameters('const', Model1, 1);
            Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,2);',data.names.X,data.names.X);
            eval(Str)
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str)
            end
            Parameters.Model1.Model = subfnSetModelParameters(Model1);
            Parameters.Model1.Outcome = data.names.M{1};
            % Model 2
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,2);',data.names.X,data.names.X);
            eval(Str)
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str)
            end
            Parameters.Model1.Model = subfnSetModelParameters(Model1);
            Parameters.Model1.Outcome = data.names.M{1};            
            
            
            % Fill in the Parameters structure with all results
            % from Model 2
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);

            Str = sprintf('Parameters.Model2.%s%d=subfnSetParameters(''%s'',Model2,1+1);',data.names.M{1},1,data.names.M{1});
            eval(Str);
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+2);',data.names.X,data.names.X);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str);
            end
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            Parameters.Model2.Outcome = data.names.M{2};
            
            % Fill in the Parameters structure with all results
            % from Model 3
            Parameters.Model3.const = subfnSetParameters('const', Model3, 1);
            Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+1);',data.names.M{2},data.names.M{2});
            eval(Str);            
            Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+2);',data.names.M{1},data.names.M{1});
            eval(Str);
            Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+3);',data.names.X,data.names.X);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+3+j);',data.names.COV{j},data.names.COV{j});
                eval(Str);
            end
            Parameters.Model3.Model = subfnSetModelParameters(Model3);
            Parameters.Model3.Outcome = data.names.Y;
            % Fill in the Parameters structure with all results
            % Model 4
            Parameters.Model4.const = subfnSetParameters('const', Model4, 1);
            Str = sprintf('Parameters.Model4.%s=subfnSetParameters(''%s'',Model4,2);',data.names.X,data.names.X);
            eval(Str)
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model4,2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str)
            end
            Parameters.Model4.Model = subfnSetModelParameters(Model4);
            Parameters.Model4.Outcome = data.names.Y;
           
            % Fill in the Parameters structure with all results
            % from the bootstrapped values

            Str = sprintf('Parameters.%s.pointEst = ParameterToBS.values(1);',ParameterToBS.names(1,:));
            eval(Str);
            Parameters.JohnsonNeyman = -99;
            Parameters.Model3.Outcome = data.names.Y;
        end
        
    case '7'
        %
        %  W  M
        %   \/ \
        %   /   \
        %  X     Y
        %
        % This is the model where the branch from X to M is moderated
        % Model 1: X to M moderated by W
        % Model 2: X and M to Y
        % Model 3: X to Y
        % Indirect effect:
        [Ndata, Nmed] = size(data.M);
        NameStruct = cell(Nmed,1);
        for j = 1:Nmed
            NameStruct{j} = sprintf('CondAB%d',j);
        end
        probeW = subfnCalculateProbeValues(data.W);
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,length(probeW)),'probeValues',zeros(1,length(probeW)),'ProbeMod',0);
        
        % Model 1
        % create the interaction terms
        Model1 = cell(Nmed);
        for j = 1:Nmed
            Model1{j} = subfnregstats(data.M(:,j),[data.X data.W data.X.*(data.W) data.COV]);
            % check each interaction term and set the flag here for whether
            % to probe the interaction.
            %if Model1{j}.tstat.pval(4) < max(data.Thresholds)
                ParameterToBS.ProbeMod = 1;
            %end
        end
        % Model 2
        % B branch model
        Model2 = subfnregress(data.Y,[data.M data.X data.COV]);
        if data.ProbeMod %m
            % check to see if the interaction is significant! This should
            % only be checked the first time through.
            % if S.tstat.pval(4) < max(data.Thresholds)
            for j = 1:Nmed
                ParameterToBS.values(j,1) = Model1{j}.beta(2)*Model2(2);
            end
            
            for k = 2:length(probeW)
                for j = 1:Nmed
                    temp = subfnregress(data.M(:,j),[data.X (data.W-probeW(k)) data.X.*(data.W - probeW(k)) data.COV]);
                    % this is the conditional parameter
                    ParameterToBS.values(j,k) = temp(2)*Model2(2);
                    ParameterToBS.probeValues(1,k) = probeW(k);
                end
            end
        else
            temp = zeros(Nmed,1);
            for j = 1:Nmed
                temp(j) = Model1{j}.beta(2)*Model2(2);
            end
            ParameterToBS.values = temp;
            ParameterToBS.probeValues = 0;
        end
        Parameters = {};
        ParameterToBS.k2 = 0;
        if PointEstFlag
            % Find the R2 increase due to the interaction
            Model2 = subfnregstats(data.Y,[data.M data.X data.COV]);
            Model3 = subfnregstats(data.Y,[data.X data.COV]);
            tcrit = tinv(1 - max(data.Thresholds)/2,length(data.X) - (4 + size(data.COV,2)));
            for j = 1:Nmed
                Parameters.Model1{j}.const = subfnSetParameters('const', Model1{j}, 1);
                Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},2);',data.names.X,data.names.X);
                eval(Str)
                Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},3);',data.names.W,data.names.W);
                eval(Str)
                Str = sprintf('Parameters.Model1{j}.%s_x_%s=subfnSetParameters(''%s_x_%s'',Model1{j},4);',data.names.X,data.names.W,data.names.X,data.names.W);
                eval(Str)
                for k = 1:size(data.COV,2)
                    Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},4+k);',data.names.COV{k},data.names.COV{k});
                    eval(Str)
                end
                Parameters.Model1{j}.Model = subfnSetModelParameters(Model1{j});
                Parameters.Model1{j}.Outcome = data.names.M{j};
                Parameters.JNvalue = subfnJohnsonNeyman(Model1{j}.beta(2),Model1{j}.covb(2,2),Model1{j}.beta(4),Model1{j}.covb(4,4),Model1{j}.covb(2,4),tcrit);
                %Parameters.JNvalue = Parameters.JNvalue.*Model2.beta(2);
            end
            
            %noIntS = subfnregstats(data.Y,[data.M data.V  data.X data.COV]);
            %diffS = subfnCalculateModelFitDiff(S,noIntS);
            %Parameters.EffOfInt = subfnSetModelParameters(diffS);
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,2);',data.names.M{j},data.names.M{j});
            eval(Str)
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,3);',data.names.X,data.names.X);
            eval(Str)
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            Parameters.Model2.Outcome = data.names.Y;
            % from Model 3
            Parameters.Model3.Outcome = data.names.Y;
            Parameters.Model3.const = subfnSetParameters('const', Model3, 1);
            Str = sprintf('Parameters.Model3.%s = subfnSetParameters(''%s'', Model3, 2);',data.names.X,data.names.X);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,3+j);',data.names.COV{j},data.names.COV{j});
                eval(Str);
                Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str);
            end
            Parameters.Model3.Model = subfnSetModelParameters(Model3);
            % calculate the Johnson-Neyman value
            %            JNvalue = subfnJohnsonNeyman(S.beta(2),S.covb(2,2),S.beta(4),S.covb(4,4),S.covb(2,4),data.tcrit);
            %             Parameters.JohnsonNeyman = JNvalue;
            %              Parameters.Bconst = subfnSetParameters('Bconst',Model2,1);
            %              Parameters.CP = subfnSetParameters('CP',Model2,Nmed+2);
            %              Parameters.C = subfnSetParameters('C',Model3,2);
            %              Parameters.Cconst = subfnSetParameters('Cconst',Model3,1);
            
            
        end
    case '14'
        %
        %     M  V
        %    / \/
        %   /   \
        %  X     Y
        %
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
        probeV = subfnCalculateProbeValues(data.V);
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,length(probeV)),'probeValues',zeros(1,length(probeV)),'ProbeMod',0);
        a = zeros(Nmed,1);
        % First, check to see if the interaction effect is significant or
        % not.
        Interaction = zeros(Ndata,Nmed);
        Model1 = cell(Nmed);
        for j = 1:Nmed
            % Use this for loop to create the interaction term for use in
            % Model2
            Interaction(:,j) = data.M(:,j).*(data.V);
            Model1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
            a(j) = Model1{j}.beta(2);
        end
        clear tempModel2
        Model2 = subfnregstats(data.Y,[data.M data.V Interaction data.X data.COV]);
        % check to see if the interaction is significant
        
        for j = 1:Nmed
        %    if Model2.tstat.pval(1+Nmed+1+j) < max(data.Thresholds)
                ParameterToBS.ProbeMod = 1;
        %    end
        end
        if data.ProbeMod
            for j = 1:Nmed
                a = Model1{j}.beta(2);
                for k = 1:length(probeV)
                    ParameterToBS.values(:,k) = a.*(Model2.beta(1+j) + Model2.beta(1+Nmed+1+1:1+Nmed+1+Nmed)*probeV(k));
                end
            end
            %                 Interaction = zeros(Ndata,Nmed);
%                 for j = 1:Nmed
%                     Interaction(:,j) = data.M(:,j).*(data.V - probeV(k));
%                     Model1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
%                     
%                 end
%                 % B branch model
%                 tempModel2 = subfnregstats(data.Y,[data.M (data.V - probeV(k)) Interaction data.X data.COV]);
%                 ParameterToBS.values(:,k) = a.*(tempModel2.beta(2:Nmed+1));
%             end
            ParameterToBS.probeValues = probeV;
        else
            ParameterToBS.values = a.*(Model2.beta(2:Nmed+1));
            ParameterToBS.probeValues = 0;
        end
        ParameterToBS.k2 = 0;
        Parameters = {};
        
        if PointEstFlag
            Interaction = zeros(Ndata,Nmed);
            Model1 = cell(Nmed);
            for j = 1:Nmed
                % Use this for loop to create the interaction term for use in
                % Model2
                Interaction(:,j) = data.M(:,j).*(data.V);
                Model1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
                a(j) = Model1{j}.beta(2);
            end
            Model2 = subfnregstats(data.Y,[data.M data.V Interaction data.X data.COV]);
            Model3 = subfnregstats(data.Y,[data.X data.COV]);
            noInt3 = subfnregstats(data.Y,[data.M data.V data.X data.COV]);
            diff3 = subfnCalculateModelFitDiff(Model3,noInt3);
            Parameters.EffOfInt = subfnSetModelParameters(diff3);
            Parameters.Model3.Model = subfnSetModelParameters(Model3);
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            % calculate the Johnson-Neyman value
            JNvalue = subfnJohnsonNeyman(Model2.beta(2),Model2.covb(2,2),Model2.beta(4),Model2.covb(4,4),Model2.covb(2,4),data.tcrit);
            Parameters.JohnsonNeyman = JNvalue;
            for j = 1:Nmed
                Parameters.Model1{j}.const = subfnSetParameters('const', Model1{j}, 1);
                Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},2);',data.names.X,data.names.X);
                eval(Str)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,j+1);',data.names.M{j},[data.names.M{j} num2str(j)]);
                eval(Str)
                Str = sprintf('Parameters.Model2.%s_x_%s=subfnSetParameters(''%s_x_%s'',Model2,Nmed+2+j);',data.names.M{j},data.names.V,data.names.M{j},data.names.V);
                eval(Str)
                for k = 1:size(data.COV,2)
                    Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},2+k);',data.names.COV{k},data.names.COV{k});
                    eval(Str)
                end
                Parameters.Model1{j}.Model = subfnSetModelParameters(Model1{j});
                Parameters.Model1{j}.Outcome = data.names.M{j};
                %                Parameters.JNvalue = subfnJohnsonNeyman(Model1{j}.beta(2),Model1{j}.covb(2,2),Model1{j}.beta(4),Model1{j}.covb(4,4),Model1{j}.covb(2,4),tcrit);
                %
                %                 Parameters.Model1.Model = subfnSetModelParameters(Model1{i});
                %                 Parameters.A{i} = subfnSetParameters('A', Model1{i},2);
                %                 Parameters.Aconst{i} = subfnSetParameters('Aconst', Model1{i},1);
                %                 Parameters.B{i} = subfnSetParameters('B', Model2,i + 1);
                %                 Parameters.Int{i} = subfnSetParameters('Int',Model2,1 + Nmed + 1 + i);
            end
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,Nmed+2);',data.names.V,data.names.V);
            eval(Str)
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+1+Nmed+1);',data.names.X,data.names.X);
            eval(Str)
            for k = 1:size(data.COV,2)
                 Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+1+Nmed+1+k);',data.names.COV{k},data.names.COV{k});
                 eval(Str)
            end
            for j = 1:Nmed
                Parameters.Model1{j}.const = subfnSetParameters('const', Model1{j}, 1);
            end
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            Parameters.Model2.Outcome = data.names.Y;
            Parameters.Model3.const = subfnSetParameters('const', Model3, 1);
            Parameters.Model3.Outcome = data.names.Y;
            Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,2);',data.names.X,data.names.X);
            eval(Str)

            for k = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,2+k);',data.names.COV{k},data.names.COV{k});
                eval(Str)
            end
%             Parameters.Bconst = subfnSetParameters('Bconst',Model2,1);
%             Parameters.V = subfnSetParameters('V',Model2,Nmed+2);
%             Parameters.C = subfnSetParameters('C',Model3,2);
%             Parameters.Cconst = subfnSetParameters('Cconst',Model3,1);
%             Parameters.CP = subfnSetParameters('Cconst',Model2,1+Nmed+1+Nmed+1);
            
        end
    case '58'
        %
        %     M
        %    / \
        %   /\ /\
        %  /  W  \
        % X       Y
        % This is a moderated mediation model. The moderator (W) affects
        % the relationship between X and M and between M and Y. The
        % conditional effect of X on Y via M is evaluated at multiple
        % moderation values. The confidence intervals for each of these
        % moderating values are calculated via bootstrapping.
        [Ndata, Nmed] = size(data.M);
        NameStruct = cell(Nmed,1);
        for j = 1:Nmed
            NameStruct{j} = sprintf('CondAB%d',j);
        end
        probeW = subfnCalculateProbeValues(data.W);
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,length(probeW)),'probeValues',zeros(1,length(probeW)),'ProbeMod',0);
        
        Model1 = cell(Nmed);
        Interaction = zeros(Ndata,Nmed);
        a = zeros(Nmed,1);
        for j = 1:Nmed
            Model1{j} = subfnregstats(data.M(:,j),[data.X data.W data.X.*(data.W) data.COV]);
            % check each interaction term and set the flag here for whether
            % to probe the interaction for BRANCH A.
            if Model1{j}.tstat.pval(4) < max(data.Thresholds)
                ParameterToBS.ProbeMod = 1;
            end
            % Use this for loop to create the interaction term for use in
            % Model2
            Interaction(:,j) = data.M(:,j).*(data.W);
            a(j) = Model1{j}.beta(2);
        end
        
        Model2 = subfnregstats(data.Y,[data.M data.W Interaction data.X data.COV]);
        % check to see if the BRANCH B interaction is significant
        for j = 1:Nmed
       %     if Model2.tstat.pval(1+Nmed+1+j) < max(data.Thresholds)
                ParameterToBS.ProbeMod = 1;
       %     end
        end
        
        % If either of the interactions are significant then probe them
         if data.ProbeMod
            % calculate the condition effect when the moderator equals zero
            %ParameterToBS.values(1,1) = (Model1)(Model2.beta(2:Nmed+1));

            

            for k = 1:length(probeW)
                Interaction = zeros(Ndata,Nmed);
                for j = 1:Nmed
                    Interaction(:,j) = data.M(:,j).*(data.W - probeW(k));
                %    tempModel1{j} = subfnregstats(data.M(:,j),[data.X data.W data.X.*(data.W - probeW(k)) data.COV]);
                    ABranchEffect(j) = Model1{j}.beta(2) + Model1{j}.beta(4).*probeW(k);
%                     BBranchEffect(j) = Model2.beta(1+j) + Model2.beta(1+Nmed+1+j).*probeW(k);
                end
                % B branch model
%                tempModel2 = subfnregstats(data.Y,[data.M (data.W-probeW(k)) Interaction data.X data.COV]);      
                for j = 1:Nmed
                     BBranchEffect(j) = Model2.beta(1+j) + Model2.beta(1+Nmed+1+j).*probeW(k);
                end
                ParameterToBS.values(:,k) = ABranchEffect.*BBranchEffect;
            end
            ParameterToBS.probeValues = probeW;
         else
            ParameterToBS.values = a.*Model2.beta(2:Nmed+1);
            ParameterToBS.probeValues = 0;
         end
         ParameterToBS.k2 = 0;
         
         if PointEstFlag
             Interaction = zeros(Ndata,Nmed);
             Model1 = cell(Nmed);
             for j = 1:Nmed
                 % Use this for loop to create the interaction term for use in
                 % Model2
                 Model1{j} = subfnregstats(data.M(:,j),[data.X data.W data.X.*(data.W) data.COV]);
                 Interaction(:,j) = data.M(:,j).*(data.W);
                 a(j) = Model1{j}.beta(2);
             end
             Model2 = subfnregstats(data.Y,[data.M data.W Interaction data.X data.COV]);
             %Model2 = subfnregstats(data.Y,[data.M data.V Interaction data.X data.COV]);
             Model3 = subfnregstats(data.Y,[data.X data.COV]);
             noInt2 = subfnregstats(data.Y,[data.M data.W data.X data.COV]);
             diff2 = subfnCalculateModelFitDiff(Model2,noInt2);
             Parameters.EffOfInt = subfnSetModelParameters(diff2);
             Parameters.Model2.Model = subfnSetModelParameters(Model2);
             Parameters.Model2.Outcome = data.names.Y;
             Parameters.Model3.Model = subfnSetModelParameters(Model3);
             Parameters.Model3.Outcome = data.names.Y;
             % calculate the Johnson-Neyman value
             JNvalue = -99;
             Parameters.JohnsonNeyman = JNvalue;
             Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1);','const','const');
             eval(Str);            
             Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+1);',data.names.W,data.names.W);
             eval(Str); 
             Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+1+Nmed+1);',data.names.X,data.names.X);
             eval(Str); 
             for j = 1:Nmed
                 Parameters.Model1{j}.const = subfnSetParameters('const', Model1{j}, 1);
                 Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},2);',data.names.X,data.names.X);%X
                 eval(Str)
                 Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},3);',data.names.W,data.names.W);%W
                 eval(Str)
                 Str = sprintf('Parameters.Model1{j}.%s_x_%s=subfnSetParameters(''%s_x_%s'',Model1{j},4);',data.names.X,data.names.W,data.names.X,data.names.W);%XxW
                 eval(Str)
                 Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,j+1);',[data.names.M{j}],[data.names.M{j}]);
                 eval(Str)
                 Str = sprintf('Parameters.Model2.%s_x_%s=subfnSetParameters(''%s_x_%s'',Model2,Nmed+2+j);',[data.names.M{j}],data.names.W,[data.names.M{j}],data.names.W);
                 eval(Str)
                 Parameters.Model1{j}.Outcome = data.names.M{j};
                 for k = 1:size(data.COV,2)
                     Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},4+k);',data.names.COV{k},data.names.COV{k});
                     eval(Str)
                 end
                 Parameters.Model1{j}.Model = subfnSetModelParameters(Model1{j});
                 Parameters.Model1{j}.Outcome = data.names.M{1};
             end
             for k = 1:size(data.COV,2)
                 Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+1+Nmed+1+k);',data.names.COV{k},data.names.COV{k});
                 eval(Str)
             end
             Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model2,1);','const','const');
             eval(Str);            
             Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,2);',data.names.X,data.names.X);
             eval(Str); 
             for k = 1:size(data.COV,2)
                 Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+k);',data.names.COV{k},data.names.COV{k});
                 eval(Str)
             end

             
%                          Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,Nmed+2);',data.names.V,data.names.V);
%             eval(Str)
%             Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+1+Nmed+1);',data.names.X,data.names.X);
%             eval(Str)
%             for k = 1:size(data.COV,2)
%                 Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model2,1+NMed+1+NMed+1+k);',data.names.COV{k},data.names.COV{k});
%                 eval(Str)
%             end
%             for j = 1:Nmed
%                 Parameters.Model1{j}.const = subfnSetParameters('const', Model1{j}, 1);
%             end
%             Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
%             Parameters.Model2.Outcome = data.names.Y;
%             Parameters.Model3.const = subfnSetParameters('const', Model3, 1);
%             Parameters.Model3.Outcome = data.names.Y;
%             Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,2);',data.names.X,data.names.X);
%             eval(Str)
% 
%             for k = 1:size(data.COV,2)
%                 Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,2+k);',data.names.COV{k},data.names.COV{k});
%                 eval(Str)
%             end

         end
         
    case '74'
         %
        %     M
        %    / \
        %   /  /\
        %  /  /  \
        % X --    Y
        %
        % This is a moderated mediation model where X is the moderator 
        % between M and Y. The conditional effect of X on Y via M is 
        % evaluated at multiple moderation values. The confidence intervals 
        % for each of these moderating values are calculated via bootstrapping.

        [Ndata Nmed] = size(data.M);
        % Check to see if the moderator, in this case X, is dicotomous.
        minX = min(data.X);
        maxX = max(data.X);
        rangeX = maxX - minX;
        if rangeX == 1
            probeX = [0 1];
        else           
            probeX = subfnCalculateProbeValues(data.X);
        end
        NameStruct = cell(Nmed,1);
        for j = 1:Nmed
            NameStruct{j} = sprintf('CondAB%d',j);
        end
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(Nmed,length(probeX)),'probeValues',zeros(1,length(probeX)),'ProbeMod',0);
        a = zeros(Nmed,1);
        % First, check to see if the interaction effect is significant or
        % not.
        Interaction = zeros(Ndata,Nmed);
        Model1 = cell(Nmed);
        for j = 1:Nmed
            % Use this for loop to create the interaction term for use in
            % Model2
            Interaction(:,j) = data.M(:,j).*(data.X);
            Model1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
            a(j) = Model1{j}.beta(2);
        end
        clear tempModel2
        tempModel2 = subfnregstats(data.Y,[data.M Interaction data.X data.COV]);
        % check to see if the interaction is significant
        
        for j = 1:Nmed
         %   if tempModel2.tstat.pval(1+Nmed+j) < max(data.Thresholds)
                ParameterToBS.ProbeMod = 1;
          %  end
        end
        if data.ProbeMod
            ParameterToBS.values(1,1) = a.*(tempModel2.beta(2:Nmed+1));

            for k = 2:length(probeX)
                Interaction = zeros(Ndata,Nmed);
                for j = 1:Nmed
                    Interaction(:,j) = data.M(:,j).*(data.X - probeX(k));
                    Model1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
                    a(j) = Model1{j}.beta(2);
                end
                % B branch model
                tempModel2 = subfnregstats(data.Y,[data.M Interaction (data.X - probeX(k)) data.COV]);
                ParameterToBS.values(:,k) = a.*(tempModel2.beta(2:Nmed+1));
            end
            ParameterToBS.probeValues = probeX;
        else
            ParameterToBS.values = a.*(tempModel2.beta(2:Nmed+1));
            ParameterToBS.probeValues = 0;
        end
        ParameterToBS.k2 = 0;
        Parameters = {};
        
        if PointEstFlag
            Interaction = zeros(Ndata,Nmed);
            Model1 = cell(Nmed);
            for j = 1:Nmed
                % Use this for loop to create the interaction term for use in
                % Model2
                Interaction(:,j) = data.M(:,j).*(data.X);
                Model1{j} = subfnregstats(data.M(:,j),[data.X data.COV]);
                a(j) = Model1{j}.beta(2);
            end
            Model2 = subfnregstats(data.Y,[data.M Interaction data.X data.COV]);
            Model3 = subfnregstats(data.Y,[data.X data.COV]);
            noInt3 = subfnregstats(data.Y,[data.M data.X data.COV]);
            diff3 = subfnCalculateModelFitDiff(Model3,noInt3);
            Parameters.EffOfInt = subfnSetModelParameters(diff3);
            Parameters.Model3.Model = subfnSetModelParameters(Model3);
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            % calculate the Johnson-Neyman value
            JNvalue = subfnJohnsonNeyman(Model2.beta(2),Model2.covb(2,2),Model2.beta(4),Model2.covb(4,4),Model2.covb(2,4),data.tcrit);
            Parameters.JohnsonNeyman = JNvalue;
            for j = 1:Nmed
                Parameters.Model1{j}.const = subfnSetParameters('const', Model1{j}, 1);
                Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},2);',data.names.X,data.names.X);
                eval(Str)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,j+1);',[data.names.M{j}],[data.names.M{j}]);
                eval(Str)
                Str = sprintf('Parameters.Model2.%s_x_%s=subfnSetParameters(''%s_x%s'',Model2,Nmed+1+j);',data.names.M{j},data.names.X,data.names.M{j},data.names.X);
                eval(Str)
                for k = 1:size(data.COV,2)
                    Str = sprintf('Parameters.Model1{j}.%s=subfnSetParameters(''%s'',Model1{j},2+k);',data.names.COV{k},data.names.COV{k});
                    eval(Str)
                end
                Parameters.Model1{j}.Model = subfnSetModelParameters(Model1{j});
                Parameters.Model1{j}.Outcome = data.names.M{1};
            end
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+Nmed+1);',data.names.X,data.names.X);
            eval(Str)
            for k = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+Nmed+Nmed+1+k);',data.names.COV{k},data.names.COV{k});
                eval(Str)
            end
            for j = 1:Nmed
                Parameters.Model1{j}.const = subfnSetParameters('const', Model1{j}, 1);
            end
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            Parameters.Model2.Outcome = data.names.Y;
            Parameters.Model3.const = subfnSetParameters('const', Model3, 1);
            Parameters.Model3.Outcome = data.names.Y;
            Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,2);',data.names.X,data.names.X);
            eval(Str)
            for k = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,2+k);',data.names.COV{k},data.names.COV{k});
                eval(Str)
            end
        end           
    case '75'
       %    M1--M2
       %   /      \
       %  /       /\
       % X--------  Y        
       %
       % Model 1: X to M1
       % Model 2: X and M1 to M2
       % Model 3: X, M1, M2 and X*M2 to Y
       % Model 4: X to Y
       
       
        [Ndata Nmed] = size(data.M);
        % Check to see if the moderator, in this case X, is dicotomous.
        minX = min(data.X);
        maxX = max(data.X);
        rangeX = maxX - minX;
        if rangeX == 1
            probeX = [0 1];
        else           
            probeX = subfnCalculateProbeValues(data.X);
        end
        
        NameStruct = {'M1M2' 'M1M2'};
        
        ParameterToBS = struct('names',char(NameStruct),'values',zeros(2,length(probeX)),'probeValues',zeros(1,length(probeX)),'ProbeMod',0);
        
        % =================================================================
        % Model 1
        tempModel1 = subfnregress(data.M(:,1),[data.X data.COV]);
        % =================================================================
        % Model 2
        tempModel2 = subfnregress(data.M(:,2),[data.M(:,1) data.X data.COV]);
        % =================================================================
        % Model 3
        Interaction = data.X.*data.M(:,2);
        tempModel3 = subfnregress(data.Y,[Interaction data.M(:,2) data.M(:,1) data.X data.COV]);
        % =================================================================
        % Model 4
        tempModel4 = subfnregress(data.Y,[data.X data.COV]);
       
        
        % =================================================================
        % Bootstrap values
        % the indirect effect which will be bootstrapped
        % X - M1 - M2 - Y
        ParameterToBS.ProbeMod = 1;
        ParameterToBS.k2 = 0;

        if data.ProbeMod
            % Where the moderating variable X has a value of zero.
            ParameterToBS.values(:,1) = tempModel1(2)*tempModel2(2)*tempModel3(2);

            for k = 2:length(probeX)
                Interaction = data.M(:,2).*(data.X - probeX(k));
                Model3 = subfnregstats(data.Y,[Interaction data.M(:,2) data.M(:,1) data.X data.COV]);
                ParameterToBS.values(:,k) = tempModel1(2)*tempModel2(2)*Model3.beta(3) + ...
                    tempModel1(2)*tempModel2(2)*Model3.beta(2)*probeX(k);
            end
            ParameterToBS.probeValues = probeX;
        else
            ParameterToBS.values(:,1) = tempModel1(2)*tempModel2(2)*tempModel3(3);
            ParameterToBS.probeValues(1) = 0;
        end
        ParameterToBS.k2 = 0;
        Parameters = {};
        
        if PointEstFlag
            Model1 = subfnregstats(data.M(:,1),[data.X data.COV]);
            Model2 = subfnregstats(data.M(:,2),[data.M(:,1) data.X data.COV]);
            Interaction = data.X.*data.M(:,2);
            Model3 = subfnregstats(data.Y,[Interaction data.M(:,2) data.M(:,1) data.X data.COV]);
            Model4 = subfnregstats(data.Y,[data.X data.COV]);
            Parameters = {};

            % Fill in the Parameters structure with all results
            % Model 1
            Parameters.Model1.const = subfnSetParameters('const', Model1, 1);
            Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,2);',data.names.X,data.names.X);
            eval(Str)
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str)
            end
            Parameters.Model1.Model = subfnSetModelParameters(Model1);
            Parameters.Model1.Outcome = data.names.M{1};
            % Model 2
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,2);',data.names.X,data.names.X);
            eval(Str)
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model1,2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str)
            end
            Parameters.Model1.Model = subfnSetModelParameters(Model1);
            Parameters.Model1.Outcome = data.names.M{1};            
            
            
            % Fill in the Parameters structure with all results
            % from Model 2
            Parameters.Model2.const = subfnSetParameters('const', Model2, 1);

            Str = sprintf('Parameters.Model2.%s%d=subfnSetParameters(''%s'',Model2,1+1);',data.names.M{1},1,data.names.M{1});
            eval(Str);
            Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+2);',data.names.X,data.names.X);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model2.%s=subfnSetParameters(''%s'',Model2,1+2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str);
            end
            Parameters.Model2.Model = subfnSetModelParameters(Model2);
            Parameters.Model2.Outcome = data.names.M{2};
            
            % Fill in the Parameters structure with all results
            % from Model 3
            Parameters.Model3.const = subfnSetParameters('const', Model3, 1);
            Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+2);',data.names.M{2},data.names.M{2});
            eval(Str);            
            Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+3);',data.names.M{1},data.names.M{1});
            eval(Str);
            Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+4);',data.names.X,data.names.X);
            eval(Str);
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model3.%s=subfnSetParameters(''%s'',Model3,1+4+j);',data.names.COV{j},data.names.COV{j});
                eval(Str);
            end
            
            Str = sprintf('Parameters.Model3.%s_x_%s=subfnSetParameters(''%s_x%s'',Model3,2);',data.names.M{2},data.names.X,data.names.M{2},data.names.X);
            eval(Str)
            
            Parameters.Model3.Model = subfnSetModelParameters(Model3);
            Parameters.Model3.Outcome = data.names.Y;
            
            % Fill in the Parameters structure with all results
            % Model 4
            Parameters.Model4.const = subfnSetParameters('const', Model4, 1);
            Str = sprintf('Parameters.Model4.%s=subfnSetParameters(''%s'',Model4,2);',data.names.X,data.names.X);
            eval(Str)
            for j = 1:size(data.COV,2)
                Str = sprintf('Parameters.Model1.%s=subfnSetParameters(''%s'',Model4,2+j);',data.names.COV{j},data.names.COV{j});
                eval(Str)
            end
            Parameters.Model4.Model = subfnSetModelParameters(Model4);
            Parameters.Model4.Outcome = data.names.Y;
           
            % Fill in the Parameters structure with all results
            % from the bootstrapped values

            Str = sprintf('Parameters.%s.pointEst = ParameterToBS.values(1);',ParameterToBS.names(1,:));
            eval(Str);
            Parameters.JohnsonNeyman = -99;
            Parameters.Model3.Outcome = data.names.Y;
        end



end
if PointEstFlag
    Parameters.names.X = data.names.X;
    Parameters.names.M{1} = data.names.M{1};
    Parameters.names.Y = data.names.Y;
    Parameters.names.V = data.names.V;
    Parameters.names.W = data.names.W;
    Parameters.ModelNum = data.ModelNum;
    Parameters.SampleSize = length(data.X);
% else
%     Parameters.Thresholds = data.Thresholds;
%     Parameters.Nboot = data.Nboot;
end
function Results = FitProcessModel(data)
% This is the program that fits the multiple regression models defined by
% this analysis.
% For the direct effects the unstandardized parameter estimates are
% calculated as beta. 
% The standardized parameter estimates are calculated as B.
% Additionally, the following are also calculated:
% t: t-statistic
% se: standard error
% p: p-value
% df: degrees of freedom for the t-test

% The interaction matrix can actually be three dimensional. This is so a
% single model effect can have multiple interaction effects in it. If there
% are multiple interactions in a full model but each path step only has one
% interaction, then the interaction matrix is only 2-D.
% This loop checks to see if there are more than on interaction effect in a
% single path step.

[N M] = size(data.data);
MaxNumberInter = 0;
if ~isempty(find(data.Inter))
    MaxNumberInter = size(data.Inter,3);
end


% Make sure the interaction matrix contains all the main effects
% if isempty(find(Inter ~= Direct.*Inter))
% The interaction terms are good

% Fit the regression equations for direct effects
Results.beta = zeros(M+1+MaxNumberInter,M);
Results.B = zeros(M+1+MaxNumberInter,M);
Results.t = zeros(M+1+MaxNumberInter,M);
Results.se = zeros(M+1+MaxNumberInter,M);
Results.p = zeros(M+1+MaxNumberInter,M);
Results.df = zeros(M+1+MaxNumberInter,M);
Results.r2 = zeros(1,M);
Results.covb = zeros(M+1+MaxNumberInter,M+1+MaxNumberInter,M);
for i = 1:M
    Col = find(data.Direct(:,i));
    if ~isempty(Col)
        % Check for interaction terms for this dependent variable
        Interaction = [];
        if ~isempty(find(data.Inter(:,i,:)))
            Interaction = zeros(N,MaxNumberInter);
            % Calculate the interaction regressor as the multiplication of
            % the regressors comprising the interaction
            for j = 1:MaxNumberInter
                InterCol = find(data.Inter(:,i,j));
                Interaction(:,j) = prod(data.data(:,InterCol),2);
            end
        end
        % Perform the regression. For maximization of speed, this
        % sub-function would need to be optimized. This could be done be
        % eliminating all the regression metrics/tests that are performed
        % and only leave the absolute minimum.
        S = ProcessRegStats(data.data(:,i),[data.data(:,Col) Interaction]);
        % Once the regression model is fit, extract the required estimated
        % values.
        %
        Results.beta([1; 1+find(data.Direct(:,i))],i) = S.beta(1:length(Col)+1);
        Results.B([1; 1+find(data.Direct(:,i))],i) = S.B(1:length(Col)+1);
        Results.t([1; 1+find(data.Direct(:,i))],i) = S.tstat.t(1:length(Col)+1);
        Results.se([1; 1+find(data.Direct(:,i))],i) = S.tstat.se(1:length(Col)+1);
        Results.p([1; 1+find(data.Direct(:,i))],i) = S.tstat.pval(1:length(Col)+1);
        % Add the covariance of the betas to the Results structure.
        % 
        % THis works. It fills in all of the covariance structure for the
        % direct effects.
        colsForParameters = [1; 1+find(data.Direct(:,i))];
        colsForInteractions = max(colsForParameters)+1:max(colsForParameters)+MaxNumberInter;
        
        for j = 1:length(colsForParameters)
            Results.covb(colsForParameters,colsForParameters(j),i) = S.covb(1:length(Col)+1,j);
        end
           
            
            
        if ~isempty(Interaction)
            % Add the interaction terms
            Results.beta(M+2:end,i) = S.beta(length(Col)+2:end);
            Results.B(M+2:end,i) = S.B(length(Col)+2:end);
            Results.t(M+2:end,i) = S.tstat.t(length(Col)+2:end);
            Results.se(M+2:end,i) = S.tstat.se(length(Col)+2:end);
            Results.p(M+2:end,i) = S.tstat.pval(length(Col)+2:end);
            % Fill in the covariance for the interaction term
            % The following does it for the interaction term and each
            % direct effect.
            for j = 1:length(colsForInteractions)
                Results.covb(colsForInteractions(j),[colsForParameters; colsForInteractions],i) = S.covb(:,length(colsForParameters)+j)';
                Results.covb([colsForParameters; colsForInteractions]',colsForInteractions(j),i) = S.covb(:,length(colsForParameters)+j);
            end
         
        end
        Results.r2(i) = S.rsquare;
    end
end


% The paths are stored in an array of cells. This is inconvenient for any
% other processes, but for now it works. The reason this was chosen was
% because the size of indirect effects all depend on the model and the
% number of paths requested. Some of the paths could have interactions
% along them. For a second order interaction this is an array of probed
% values stored in the path. If the interaction is third order then the
% probed path becomes a plane. Therefore, to accomodate this variability in
% the number of dimensions an array of cells was chosen. 
NumberOfPaths = size(data.Paths,3);
Results.Paths = cell(NumberOfPaths,1);
Results.PathsSE = cell(NumberOfPaths,1);
Results.PathsTnorm = cell(NumberOfPaths,1);
Results.ProbeValues = cell(MaxNumberInter,MaxNumberInter);
[OutPath, OutSE, ProbeValues] = subfnCalculatePathSE(Results, data);
Results.Paths{1} = OutPath;
Results.PathsSE{1} = OutSE;
Results.PathsTnorm{1} = OutPath./OutSE;
Results.ProbeValues{1} = ProbeValues;
%[PathsParameters, PathsStandardErrors] = subfnCalculatePathSE(Results, data)

%% add a row of zeros so that the indexing matchs with the beta matrix
% Paths = [zeros(1,M,NumberOfPaths); data.Paths];
% 
% for j = 1:NumberOfPaths
%     % step through the path
%     ResultPath = 1;
%     ResultPathbetas = zeros(max(max(Paths(:,:,j))),1);
%     ResultPathSE = zeros(max(max(Paths(:,:,j))),1);
%     
%     for i = 1:max(max(Paths(:,:,j)))
%         
%         index = find(Paths(:,:,j)==i);
%         [Row Col] = ind2sub(size(Paths),index);
%         % is there an interaction on this step?
%         tempInter = data.Inter(:,Col);
%         
%         if isempty(find(tempInter))
%             ResultPath = ResultPath.*Results.beta(Row,Col);
%             ResultPathbetas(i) = Results.beta(Row,Col);
%             ResultPathSE(i) = Results.se(Row,Col);
%         else
%             % There is an interaction here
%             InteractionComponent = Results.beta(Row,Col);
%             % which variable do you probe?
%             
%             % find the probe values,this even works for higher
%             % order interactions, I think.
%             F = find(tempInter) + 1;
%             % probe the one NOT in the path
%             % Check to see if the interaction effect is even in the path
%             if ~isempty(find(F == Col))
%                 F(find(F == Col)) = 0;
%                 % do not probe a
%                 Fmod = F(find(F)) - 1;
%                 Moderators = data.data(:,Fmod);
%                 if length(unique(Moderators)) == 2
%                     probeMod = unique(Moderators);
%                     probeValues = probeMod;
%                 else
%                     %probeMod = prctile(Moderators,[10:10:90]);
%                     probeMod = prctile(Moderators,[5:5:95]);
%                     probeValues = [zeros(size(probeMod,1),1) probeMod];
%                     %  probeValues = probeValues(:,1)*probeValues(:,2)';
%                 end
%                 InteractionComponent = InteractionComponent + Results.beta(M+1+1,Col).*probeValues;
%                 ResultPath = ResultPath.*InteractionComponent;
%             else
%                 ResultPath = ResultPath.*Results.beta(Row,Col);
%             end
%         end
%     end
%     % Store the probe values and the resultant path values.
%     Results.ProbeValues{j} = probeValues;
%     Results.Paths{j} = ResultPath;
%     if length(ResultPathSE) == 2
%         Results.PathsSE{j} = sqrt(sum((ResultPathbetas.^2).*flipud((ResultPathSE.^2))));
%     elseif length(ResultPathSE) == 3
%         Results.PathsSE{j} = sqrt((prod(ResultPathbetas(1:2).^2).*(ResultPathSE(3).^2)) + ...
%         (prod(ResultPathbetas([1 3]).^2)).*(ResultPathSE(2).^2) + ...
%         (prod(ResultPathbetas([2 3]).^2)).*(ResultPathSE(1).^2));
%     else
%         Results.PathsSE{j} = 0;
%     end
%     Results.PathsTnorm{j} = ResultPath/Results.PathsSE{j};
% end
% 
% 

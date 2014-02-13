function Results = WIPsubfnFitModel(data)
[N M] = size(data.data);
MaxNumberInter = 0;
if ~isempty(find(data.Inter))
    MaxNumberInter = size(data.Inter,3);
end

% Make sure the interaction matrix contains all the main effects
%if isempty(find(Inter ~= Direct.*Inter))
% The interaction terms are good
%end
% Fit the regression equations for direct effects
Results.beta = zeros(M+1+MaxNumberInter,M);
Results.B = zeros(M+1+MaxNumberInter,M);
Results.t = zeros(M+1+MaxNumberInter,M);
Results.se = zeros(M+1+MaxNumberInter,M);
Results.p = zeros(M+1+MaxNumberInter,M);
Results.df = zeros(M+1+MaxNumberInter,M);
for i = 1:M
    Col = find(data.Direct(:,i));
    if ~isempty(Col)
        % Check for interaction terms for this dependent variable
        Interaction = [];
        if ~isempty(find(data.Inter(:,i,:)))
            Interaction = zeros(N,MaxNumberInter);
            for j = 1:MaxNumberInter
                InterCol = find(data.Inter(:,i,j));
                Interaction(:,j) = prod(data.data(:,InterCol),2);
            end
        end
        % b = subfnregress(Data(:,i),[Data(:,Col) Interaction]);
        S = subfnregstats(data.data(:,i),[data.data(:,Col) Interaction]);
        Results.beta([1; 1+find(data.Direct(:,i))],i) = S.beta(1:length(Col)+1);
        Results.B([1; 1+find(data.Direct(:,i))],i) = S.beta(1:length(Col)+1);
        Results.t([1; 1+find(data.Direct(:,i))],i) = S.tstat.t(1:length(Col)+1);
        Results.se([1; 1+find(data.Direct(:,i))],i) = S.tstat.se(1:length(Col)+1);
        Results.p([1; 1+find(data.Direct(:,i))],i) = S.tstat.pval(1:length(Col)+1);
        % df([1; 1+find(Direct(:,i))],i) = S.tstat.dfe;
        if ~isempty(Interaction)
            % Add the interaction terms
            Results.beta(M+2:end,i) = S.beta(length(Col)+2:end);
            %B(M+2:end,i) = B(length(Col)+2:end);
            Results.t(M+2:end,i) = S.tstat.t(length(Col)+2:end);
            Results.se(M+2:end,i) = S.tstat.se(length(Col)+2:end);
            Results.p(M+2:end,i) = S.tstat.pval(length(Col)+2:end);
            %    df(M+2:end,i) = S.tstat.dfe;
            
        end
    end
end



NumberOfPaths = size(data.Paths,3);
Results.IndirectEffect = cell(size(data.Paths,3),1);
% add a row of zeros so that the indexing matchs with the beta matrix
Paths = [zeros(1,M,NumberOfPaths); data.Paths];
for j = 1:NumberOfPaths
    % cycle over the steps in this path
    Results.IndirectEffect{j} = 1;
    IndirectEffect2D = 1;
    for i = 1:max(max(data.Paths(:,:,j)))
        
        index = find(Paths(:,:,j)==i);
        [Row Col] = ind2sub(size(Paths),index);
        
        % are there any interactions at this point in the path?
        % cycle over all interactions in model
        for k = 1:MaxNumberInter
            
            if ~isempty(find(data.Inter(:,Col,k)))
                % yes there is an interaction
                % probe the interaction
                
                % find the OTHER terms that interact with the path of
                % interest
                tempInter = data.Inter(:,Col,k);
                tempInter(Row) = 0;
                % find the probe values,this even works for higher
                % order interactions, I think
                ThisStepInter = Results.beta(Row,Col) + Results.beta(M+1+k,Col).*[0; prod(prctile(data.data(:,find(tempInter)),[10:10:90]),2)];
                if length(IndirectEffect2D) == 1
                    IndirectEffect2D=ThisStepInter;
                else
                    IndirectEffect2D = IndirectEffect2D*ThisStepInter';
                end
                Results.IndirectEffect{j} = Results.IndirectEffect{j}.*ThisStepInter;
            else
                % no interaction at this step
                Results.IndirectEffect{j} = Results.IndirectEffect{j}.*Results.beta(Row,Col);
                IndirectEffect2D = IndirectEffect2D.*Results.beta(Row,Col);
            end
            
            
        end
        
    end
end

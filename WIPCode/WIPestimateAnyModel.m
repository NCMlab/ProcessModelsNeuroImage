clear
N = 100;
M = 3;
Data = randn(N,M);
% Create the direct effects model
Direct = zeros(M,M);
Direct([1],2) = 1;
Direct([1 2],3) = 1;

% Interactions
Inter = zeros(M,M);
%Inter([4 5],3,1) = 1;
Inter([1 2 ],2,1) = 1;

% Estimate the paths
Paths = zeros(M,M,1);
Paths(1,2,1) = 1;
Paths(2,3,1) = 2;


Direct
Inter
Paths
%
MaxNumberInter = 0;
if ~isempty(find(Inter))
    MaxNumberInter = size(Inter,3);
end

% Make sure the interaction matrix contains all the main effects
%if isempty(find(Inter ~= Direct.*Inter))
% The interaction terms are good
%end
% Fit the regression equations for direct effects
beta = zeros(M+1+MaxNumberInter,M);
B = zeros(M+1+MaxNumberInter,M);
t = zeros(M+1+MaxNumberInter,M);
se = zeros(M+1+MaxNumberInter,M);
p = zeros(M+1+MaxNumberInter,M);
df = zeros(M+1+MaxNumberInter,M);
for i = 1:M
    Col = find(Direct(:,i));
    if ~isempty(Col)
        % Check for interaction terms for this dependent variable
        Interaction = [];
        if ~isempty(find(Inter(:,i,:)))
            Interaction = zeros(N,MaxNumberInter);
            for j = 1:MaxNumberInter
                InterCol = find(Inter(:,i,j));
                Interaction(:,j) = prod(Data(:,InterCol),2);
            end
        end
       % b = subfnregress(Data(:,i),[Data(:,Col) Interaction]);
        S = subfnregstats(Data(:,i),[Data(:,Col) Interaction]);
        beta([1; 1+find(Direct(:,i))],i) = S.beta(1:length(Col)+1);
        B([1; 1+find(Direct(:,i))],i) = S.beta(1:length(Col)+1);
        t([1; 1+find(Direct(:,i))],i) = S.tstat.t(1:length(Col)+1);
        se([1; 1+find(Direct(:,i))],i) = S.tstat.se(1:length(Col)+1);
        p([1; 1+find(Direct(:,i))],i) = S.tstat.pval(1:length(Col)+1);
       % df([1; 1+find(Direct(:,i))],i) = S.tstat.dfe;
        if ~isempty(Interaction)
            % Add the interaction terms
            beta(M+2:end,i) = S.beta(length(Col)+2:end);
            %B(M+2:end,i) = B(length(Col)+2:end);
            t(M+2:end,i) = S.tstat.t(length(Col)+2:end);
            se(M+2:end,i) = S.tstat.se(length(Col)+2:end);
            p(M+2:end,i) = S.tstat.pval(length(Col)+2:end);
        %    df(M+2:end,i) = S.tstat.dfe;
            
        end
    end
end
beta


NumberOfPaths = size(Paths,3);
IndirectEffect = cell(size(Paths,3),1);
% add a row of zeros so that the indexing matchs with the beta matrix
Paths = [zeros(1,M,NumberOfPaths); Paths];
for j = 1:NumberOfPaths
    % are there any interactions along this path?
    if isempty(find(Paths(2:end,:,j).*Inter))
        % no interactions
        IndirectEffect{j} = 1;
        for i = 1:max(max(Paths(:,:,j)))
            index = find(Paths(:,:,j)==i);
            [Row Col] = ind2sub([M+1,M],index);
            IndirectEffect{j} = IndirectEffect{j}.*beta(index);
        end
    else
        IndirectEffect{j} = ones(10,1);
        for i = 1:max(max(Paths(:,:,j)))
            index = find(Paths(:,:,j)==i);
            [Row Col] = ind2sub([M+1,M],index);
            % is there an interaction term for this column?
            % Only probe an interaction if there is an interaction at some
            % point in the path
            if MaxNumberInter
                for k = 1:MaxNumberInter
                    if ~isempty(Inter(Row,Col,k))
                        % find the OTHER terms that interact with the path of
                        % interest
                        tempInter = Inter(:,Col,k);
                        tempInter(Row) = 0;
                        % find the probe values,this even works for higher
                        % order interactions, I think
                        IndirectEffect{j} = IndirectEffect{j}.*beta(index) + beta(M+1+k,Col).*[0; prod(prctile(Data(:,find(tempInter)),[10:10:90]),1)'];
                    end
                end
            
            end
            
        end
    end
end
IndirectEffect
%%
data = []
data.names.X = 'X';
data.names.M{1} = 'M';
data.names.Y = 'Y';
data.names.V='';
data.names.W='';
data.X = Data(:,1);
data.M = Data(:,2);
data.Y = Data(:,3);
data.COV = [];
data.V = [];
data.ModelNum = '74'
data.ProbeMod=1
PointEstFlag = 1;
[ParameterToBS Parameters] = subfnProcessModelFit(data,PointEstFlag)


Parameters.Model1{1}



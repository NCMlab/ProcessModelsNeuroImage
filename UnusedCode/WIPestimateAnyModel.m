clear
N = 100;
M = 5;
data.data = randn(N,M);
data.data(:,3) = data.data(:,3) + data.data(:,2).*0.1;
data.STRAT = round(rand(N,1));
% Create the direct effects model
Direct = zeros(M,M);
Direct([1 4 5],2) = 1;
Direct([1 2 4 5],3) = 1;

% Interactions
Inter = zeros(M,M);

%Inter([4 5],3,1) = 1;
%Inter([1 4],2,1) = 1;
%Inter([2 4],3,1) = 1;
%Inter([1 2 4],3,1) = 1;

% Estimate the paths
Paths = zeros(M,M,1);
Paths(1,2,1) = 1;
Paths(2,3,1) = 2;


data.Direct = Direct;
data.Inter = Inter;
data.Paths = Paths;
Direct
Inter
Paths
%%

Results = WIPsubfnFitModel(data)
Results.beta
Results.Paths
%%
Nboot = 50;
tic
BootStrap = WIPsubfnBootStrap(data,Nboot);
toc
%%


Results.beta
Results.B
Results.t


X = randn(N,1);
Y = randn(N,1) + 0.5.*X;
S = subfnregstats(Y,X)





    %%
    data = []
    data.names.X = 'X';
    data.names.M{1} = 'M';
    data.names.Y = 'Y';
    data.names.V='';
    data.names.W='';
    data.X = data(:,1);
    data.M = data(:,2);
    data.Y = data(:,3);
    data.COV = [];
    data.V = [];
    data.ModelNum = '74'
    data.ProbeMod=1
    PointEstFlag = 1;
    [ParameterToBS Parameters] = subfnProcessModelFit(data,PointEstFlag)
    
    
    Parameters.Model1{1}
    
    

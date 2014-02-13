clear
N = 100;
M = 4;
data.data = randn(N,M);
data.STRAT = round(rand(N,1));
% Create the direct effects model
Direct = zeros(M,M);
Direct([1 4],2) = 1;
Direct([1 2 4],3) = 1;

% Interactions
Inter = zeros(M,M);
%Inter([4 5],3,1) = 1;
%Inter([1 4],2,1) = 1;
Inter([2 4],3,1) = 1;
%Inter([1 4],2,1) = 1;

% Estimate the paths
Paths = zeros(M,M,1);
Paths(1,2,1) = 1;
Paths(2,3,1) = 2;


data.Direct = Direct;
data.Inter = Inter;
data.Paths = Paths;
%

Results = WIPsubfnFitModel(data)

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
    
    

clear
N = 113;
Nvar = 3;
RawData = [round(rand(N,1)) randn(N,Nvar-1)];
%RawData = randn(N,Nvar);
% prepare data

% prepare names
names = cell(1,Nvar);
names{1} = 'A';
names{2} = 'B';
names{3} = 'C';
names{4} = 'D';
% prepare model effects
% Direct
Direct = zeros(Nvar);
Direct(1,[2 3]) = 1;
Direct(2,3) = 1;
% Direct(1,[2 3 4]) = 1;
% Direct(2,[3 4]) = 1;
% Direct(3,[4]) = 1;

% interactions
Inter = zeros(Nvar);
%Inter([1 3],4) = 1;
Inter([1 2],3) = 1;
% Path(s) to test
Paths = zeros(Nvar);
Paths(1,2) = 1;
Paths(2,3) = 2;
%Paths(3,4) = 3;

%Paths(1,2,2) = 2;
%Paths(2,3,2) = 1;

Model = {};
Model.data = RawData;
Model.names = names;
Model.Direct = Direct;
Model.Inter = Inter;
Model.Paths = Paths;
 
Model.N = N;
Model.Nvar = Nvar;
Model.Nboot = 500;
Model.Stratify = [];
Model.Nperm = 0;
Model.NJobSplit = 1;
Model.Thresh = [0.05];

% View the model using the biograph toolbox
% B = biograph(Model.Direct,Model.names)
% view(B);
% Use annotation


%
Results = WIPsubfnProcessData(Model)
WIPPrintResults(Model,Results)
%%
% How many interaction effects are there?
D = zeros(Nvar*2,Nvar*2);
D(1,[4 6]) = 1;
D(2,5) = 1;
D(4,2) = 1;
D(5,3) = 1;
D(6,3) = 1;
names = {'A' 'B' 'C' 'o1' 'o2' 'o3'}
B = biograph(D,names);
view(B)

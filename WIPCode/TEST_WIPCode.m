N = 100;
Nvar = 3;
RawData = randn(N,Nvar);
% prepare data

% prepare names
names = cell(1,Nvar);
names{1} = 'A';
names{2} = 'B';
names{3} = 'C';

% prepare model effects
% Direct
Direct = zeros(Nvar);
Direct(1,[2 3]) = 1;
Direct(2,[3]) = 1;
% interactions
Inter = zeros(Nvar);
% Path(s) to test
Paths = zeros(Nvar);
Paths(1,2) = 1;
Paths(2,3) = 2;


Model = {};
Model.data = RawData;
Model.names = names;
Model.Direct = Direct;
Model.Inter = Inter;
Model.Paths = Paths;
 
Model.N = N;
Model.Nvar = Nvar;
Model.Nboot = 5000;
Model.Stratify = [];
Model.Nperm = 0;
Model.NJobSplit = 1;
Model.Thresh = [0.025 0.005];


Results = WIPsubfnProcessData(Model)
WIPPrintResults(Model,Results)
function TestModel74

load Model74TestData
NSub = length(Sex);
MissingValues = ones(NSub,1);
MissingValues(find(Memory == -9999)) = 0;

Include = find(MissingValues);

data = {};
data.names = {};
data.names.X = 'AgeGroup';
data.names.M = {'RightHipp'};
data.names.Y = 'Memory';
data.names.V = 'V';
data.names.W = 'W';
data.names.COV = {'sex' 'TIV' 'GMV'}
data.STRAT = [];
data.COV = [];%randn(N,2);
data.V = [];
data.W = [];
data.Q = [];
data.R = [];
data.ModelNum = '74';
data.Thresholds = [0.05];
data.Indices = 1;
data.Nboot = 2000;

% Model of interest

data.X = AgeGroup(Include);
data.Y = Memory(Include);
data.M = RHippoVol(Include);
data.COV = [Sex(Include) TIV(Include) TotalGray(Include)]


tic
Parameters = subfnVoxelWiseProcessBatch(data);
toc


subfnPrintResults(Parameters{1})

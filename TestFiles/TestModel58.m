clear
load('/share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/PracticeData/ModMeddata.mat')



data = {};
data.names = {};
data.names.X = 'X';
data.names.M = {'M'};
data.names.Y = 'Y';
data.names.V = 'V';
data.names.W = 'V';
data.names.Q = '';
data.names.R = '';
data.names.COV = {'cov1' 'cov2'};

data.STRAT = [];
data.COV = COV;
data.W = V;
data.Q = [];
data.R = [];

data.ModelNum = '14';
data.Thresholds = [0.0001];
data.Indices = 1;
data.Nboot = 2000;

% Model of interest
data.X = X;
data.Y = Y;
data.M = M;
data.V = V;
% Add the bootstrap resamples to the input data
NSub = N;

%data.Resamples = subfuncBootStrapResamples(data.Nboot,NSub,data.STRAT);

Parameters = subfnVoxelWiseProcessBatch(data);
subfnPrintResults(Parameters{1})
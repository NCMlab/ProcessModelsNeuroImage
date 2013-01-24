clear
load('/share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/PracticeData/ModMeddata.mat')



data = {};
data.names = {};
data.names.X = 'A';
data.names.M = {'B'};
data.names.Y = 'C';
data.names.V = 'V';
data.names.W = 'W';
data.names.Q = '';
data.names.R = '';
data.names.COV = {'cov1' 'cov2'};

data.STRAT = [];
data.COV = COV;
data.W = V;
data.Q = [];
data.R = [];

data.ModelNum = '58';
data.Thresholds = [0.05];
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
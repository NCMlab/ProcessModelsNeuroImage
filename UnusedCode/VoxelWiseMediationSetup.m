% Need PInd which is the index of all voxels to include in the analysis
% This is taken from a mask image
%
% Need data in arrays: [Nsub x Nmed x Nvoxels] where the Nvoxels is equal
% to the length of the voxels in the mask image.
%
AllData = {};
% fill in names
AllData.names = {};
AllData.names.Y = 'CognitiveMeasure';
AllData.names.X = 'AgeGroup';
AllData.names.M = {'BrainMeasure'}; % Notice that this is a cell array which allows for multiple mediators
AllData.names.V = '';
AllData.names.W = '';
AllData.names.Q = '';
AllData.names.R = '';
% fill in data
AllData.Y = sRT;
AllData.M = GMdata;
AllData.X = Group;
AllData.COV = [];
AllData.STRAT = []; % stratify the resmpling? This is a good idea if you do not have the same N in each group.
AllData.V = [];
AllData.W = [];
AllData.Q = [];
AllData.R = [];
AllData.Indices = PInd;
AllData.DIM = M.dim;
Nboot = 5000;
ModelNum = '4';
Thresholds = [0.05 0.01 0.005];
NJobSplit = 150;
[Nsub Nmed Nvoxels] = size(AllData.M);

% Create a parameter file
AnalysisParameters = {};
AnalysisParameters.BaseDir = '/share/data/users/js2746_Jason/CogReserveAnalyses';
AnalysisParameters.Nsub = Nsub;
AnalysisParameters.Nmed = Nmed;
AnalysisParameters.Nvoxels = Nvoxels;
AnalysisParameters.NJobSplit = NJobSplit;
AnalysisParameters.Nboot = Nboot;
AnalysisParameters.ModelNum = ModelNum;
AnalysisParameters.Thresholds = Thresholds;
AnalysisParameters.Indices = AllData.Indices;
AnalysisParameters.names = AllData.names;
AnalysisParameters.Tag = ['X_M_Y'];

subfnRunVoxelwiseProcess(AllData,AnalysisParameters);

% after the above finishes, then run this:
WriteOutResults;
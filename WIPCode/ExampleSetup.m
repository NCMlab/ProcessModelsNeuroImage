addpath /share/studies/CogRes/Scripts/ProcessModelsNeuroImage/FinalCode
clear
BaseDir = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/SampleSize140';
cd(BaseDir)
load('Age_GM_ECF_PERF_LargerSample');

Nsub = length(AgeGroup1old);

Pmask = fullfile(BaseDir,'mask','maskIncGrey50.nii');
Vmask = spm_vol(Pmask);
Imask = spm_read_vols(Vmask);
Indices = find(Imask);
Nvoxels = length(Indices);

% Load up Functional Data
fprintf(1,'Loading fMRI data...\n');
FMRI = zeros(Nsub,1,Nvoxels);
for i = 1:Nsub
    CurrentPath = FMRIPATHS{i};
    CurrentPath = strrep(CurrentPath,'spmT_','con_');
    
    V = spm_vol(CurrentPath);
    I = spm_read_vols(V);
    FMRI(i,1,:) = I(Indices);
end

%% The aim of this example is to demonstrate a setup of this software.
% This will be an incremental step example. 
% First, will be a simple mediation analysis where the effect of age group
% on a behavioral performance measure is mediated by a voxel-wise measure
% of grey matter density.

% To begin, start by loading up the brain data. 
% This can be done by selecting the images using the spm_image function for
% file selection.  Or as below by loading up a structure of file paths.
% Determine the numner of subjects selected.
% NSUB

% A mask image is needed. From this mask image the indices of voxels to
% include is determined.
% MASK

% Load up Structural Data
fprintf(1,'Loading structural data...\n');
% allocate memory for storing the image data
STRUCTURE = zeros(Nsub,1,Nvoxels);
% cycle over the number of subjects
for i = 1:Nsub
    % read each header
    V = spm_vol(STRUCTUREPATHS{i});
    % read each image
    I = spm_read_vols(V);
    
    STRUCTURE(i,1,:) = I(Indices);
end

% For this example the behavioral data is a single number per person.
% Behavioral Data
BEHAVIOR = ECFDualCORmedRT - ECFSingCORmedRT;


fprintf(1,'Done prepapring data.\n');

%% General Setup
% Create a structure where each cell is a node in the path diagram. This
% could be a voxel-wise matrix or a vector.
data = {};
data{1} = AgeGroup1old;
data{2} = squeeze(STRUCTURE);
data{3} = squeeze(FMRI);
data{4} = BEHAVIOR;
data{5} = nWBV;
data{6} = Sex;
% Create a structure of corresponding names for each cell in the data
% structure
Names = {};
Names{1} = 'AgeGr';
Names{2} = 'pGM';
Names{3} = 'fMRI';
Names{4} = 'medRTSC';
Names{5} = 'nWBV';
Names{6} = 'Sex';

% How many variables were entered
Nvar = length(data);
% Right now the statistics can be performed at the voxel-level using bias
% correctsed accelerated confidence intervals determined from bootstrap
% resampling.
% Are there any voxel-wise bootstrap resamplings?
Nboot = 0;
% An alternative to voxel-wise statistics is a map-wise statistics based
% off of the maximum statistic approach from a series of permutation
% resamples. This approach performs the voxel-wise calculations, then
% identifiies the value of the maximal statistic in the image. This is
% redone with each permutation to create a distribution of maximal
% statistical values. These values are sorted and the percentiles are
% determined based onthe number of permutations and the thresholds. Note,
% two-tailed thresholds are used.
%
% Are there any permutaion resamples to perform?
Nperm = 5000;
% The job split variable is how many jobs this analysis is split into for
% sending to a comuputer cluster environment. Note that the more splits
% does not mean faster computing because there is a limitationon the number
% of available cluster nodes. When the number of job splits exceeds the
% number of available jobs then jobs are placed inthe queue and need to
% wait.
NJobSplit = 500;
% Specify the thresholds used in the analysis. 
% Some of the images and results in the mediation analysis perform
% parametric tests where the threshold can be more dynamically changed.
% However, for  estimation of the significance of the paths resampling
% methods are required and the probabilities need to be calculated at the
% time of the bootstrapping and therefore, made a prior.
Thresh = [0.025 0.005];
% Create a structure which will contain all information for this analysis.
ModelInfo = {};
ModelInfo.Names = Names;
ModelInfo.data = data;
ModelInfo.Nboot = Nboot;
ModelInfo.Nperm = Nperm;
ModelInfo.BaseDir = BaseDir;

ModelInfo.Indices = Indices;
ModelInfo.NJobSplit = NJobSplit;
ModelInfo.Thresholds = Thresh;
% Startification is used when the resamples are created and needs to be a
% binomial parameter for right now. The use of a stratification variable is
% so that when the reampling is performed each resample maintains the
% number of subjects as in the stratification parameter. This is most
% applicable when there are multiple groups with different sample sizes.
ModelInfo.STRAT = [];
ModelInfo.NSub = size(data{1},1);
ModelInfo.Nvar = Nvar;

%% Model Specific
% Model 75
Model75 = ModelInfo;
Model75.Tag = 'Model75_WithCov';
% The modeling is specified using a series of matrices. The first is the
% DIRECT matrix. This matrix is square with dimension based onthe number of
% variables in the model. The rows and columns refer to the variabvles in
% the model in the same order as they are specified in the data matrix.
% Since the path analyses are in fact a series of regression models these
% models get specified in DIRECT matrix. Interaction or moderating effects
% get specified in the INTER(actions) matrix below.

% A simple mediation model would be: 
% VAR_2 = VAR_1
% VAR_3 = VAR_1 + VAR_2

Direct = zeros(Nvar);
Direct(1,[2 3]) = 1;
Direct(2,[3]) = 1;
% To specify interactions in the model create another matrix the same size
% as the DIRECT matrix.
% An interaction between two variables in predicting a third, the rows
% of the interacting variables will both have a value of one within the
% column of the variable they are predicting. BE sure to include in the
% DIRECT matrix the main effects for each of these variables on the
% predicted variable.
Inter = zeros(Nvar);
% The paths to be tested are included in another matrix the same size of
% the DIRECT matrix. One difference here is that the PATHS matrix may have
% a third dimension to test multiple paths within a single model. The steps
% along a path are specified with integers, i.e. step 1, step 2 ...
Paths = zeros(Nvar);
Paths(1,2) = 1;
Paths(2,3) = 2;
Paths(3,4) = 3;

Model75.Direct = Direct;
Model75.Inter = Inter;
Model75.Paths = Paths;

WIPsubfnPrepareDataForProcessPERMUTE(Model75)

% Model 6
Model6 = ModelInfo;
Model6.Tag = 'Model6_WithCov';
Direct = zeros(Nvar);
Direct(1,[2 3 4]) = 1;
Direct(2,[3 4]) = 1;
Direct(3,[4]) = 1;
Direct([5 6],[2 3 4]) = 1;
Inter = zeros(Nvar);
%Inter([1 3],4) = 1;
Paths = zeros(Nvar);
Paths(1,2) = 1;
Paths(2,3) = 2;
Paths(3,4) = 3;

Model6.Direct = Direct;
Model6.Inter = Inter;
Model6.Paths = Paths;
WIPsubfnPrepareDataForProcessPERMUTE(Model6)

% TODO 
% The program that prepares the data to be submitted to the cluster also
% needs to decide whether to launch the permuattion test split or the
% bootstrap split. 
% In the case of the permutation test a single copy of the data is saved.
% For bootstrap the full data set is saved and then a broken up copy of
% data is also saved. Each of these is a single "chunk" of the data based
% onthe number of job splits. There is the option of not creating and
% saving all these samll chunks but to pass the path to the full data set
% and only break up the indices to analyze by each compute node. The only
% downside I see to this is the potential to use a lot of RAM memory buy
% loading up the very large data set for each compute node. And since each
% computer in the cluster has muliple nodes, the full data set will be
% loaded up by each node on each computer. It is possible to load the full
% data set and then exract the data chunk of interest and rgen remove the
% non-included data from memory. The problem is that there will be
% transient times of large amounts of data in memory. I am afraid to
% overload a single computer even transiently.
% 
% 
% 
% WIPsubfnPrepareDataForProcessPERMUTE
% % THis performs the permutation testing approach
% WIPsubfnVoxelWiseProcessPERMUTE(InDataPath,count,Nperm)
% % This applies the BCaCI to each voxel
% Results = WIPsubfnProcessData(Model)
%     
%     BootStrap = WIPsubfnBootStrap(Model,Model.Nboot,FieldNames);
% 
%     % Perform the jack-knife step
%     JackKnife = WIPJackKnife(Model,Results,FieldNames);
%     
%     % Calculate the BCaci values for each parameter
%     Results.BCaCI = WIPsubfnCreateBCaCI(Results,BootStrap,JackKnife,Model.Thresh);
% 
% % THis does the actual regression model fitting.
% Results = WIPsubfnFitModel(data)
% 

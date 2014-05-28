%% The aim of this example is to demonstrate a setup of this software.
% This will be an incremental step example. 
% First, will be a simple mediation analysis where the effect of age group
% on a behavioral performance measure is mediated by a voxel-wise measure
% of grey matter density.

% To begin, start by loading up the brain data. 
P = spm_select(Inf,'image','Select the imaging data which will serve as the mediator.');
V = spm_vol(P);
I = spm_read_vols(V);
% The format of the input data is needed for later writing of the output
% images
DataHeader = V(1);

% A mask image is needed. From this mask image the indices of voxels to
% include is determined.
% MASK
MaskFlag = 0;
while MaskFlag == 0
    Pmask = spm_select(1,'image','Select the mask image');
    Vmask = spm_vol(Pmask);
    Imask = spm_read_vols(Vmask);
    % make sure the selected mask has the same dimensions as the data
    if false(Vmask.dim == V(1).dim)
        errordlg('The mask is not the same size as the selected data');
    end
    % make sure the mask image is really a mask image
    if length(unique(Imask)) > 2
        errordlg('The mask image selected has more then two unique values.');
    else
        % Find the indices of the voxels included inthe mask
        Indices = find(Imask);
        clear Imask Vmask
        MaskFlag = 1;
    end
end
% determine the size of the data
Nvoxels = length(Indices);
Nsub = size(I,4);
% restructure the imaging data to be a Nvoxels X NSub matrix
NIData = zeros(Nsub, Nvoxels);
for i = 1:Nsub
    temp = I(:,:,:,i);
    NIData(i,:) = temp(Indices)';
end
% Clean up the memory
clear temp I 

% For this example the behavioral data is a single number per person.
% Behavioral Data
BEHAVIOR = randn(Nsub,1);
AgeGroup = round(rand(Nsub,1));

COVARIATES = randn(Nsub,2);

fprintf(1,'Done preparing data.\n');

%% General Setup
% The base directory is the folder containg the folder of results. The
% results for a specific analysis will all be in folders contained within
% this base directory. The output folders are named according to a
% user specified "Tag" name.

% Check to see whether this example is running on my laptop or desktop.
% This will have to be changed for any specific user.
if ismac
    BaseDir = '/Users/jason/Dropbox/SteffenerColumbia/Projects/TestProcessCode';
elseif ispc
    BaseDir = 'C:\Users\js2746\Dropbox\SteffenerColumbia\Projects\TestProcessCode';
elseif isunix
    BaseDir = '/home/js2746/DropBox/SteffenerColumbia/Projects/TestProcessCode';
end

% Create a structure where each cell is a node in the path diagram. This
% could be a voxel-wise matrix or a vector.
data = {};
data{1} = AgeGroup;
data{2} = NIData;
data{3} = BEHAVIOR;

% Create a structure of corresponding names for each cell in the data
% structure
Names = {};
Names{1} = 'AgeGr';
Names{2} = 'Brain';
Names{3} = 'Behavior';

% How many variables were entered
Nvar = length(data);

% Right now the statistics can be performed at the voxel-level using bias
% correctsed accelerated confidence intervals determined from bootstrap
% resampling.
% Are there any voxel-wise bootstrap resamplings?
Nboot = 100;

% An alternative to voxel-wise statistics is a map-wise statistics based
% off of the maximum statistic approach from a series of permutation
% resamples. This approach performs the voxel-wise calculations, then
% identifiies the value of the maximal statistic in the image. This is
% redone with each permutation to create a distribution of maximal
% statistical values. These values are sorted and the percentiles are
% determined based onthe number of permutations and the thresholds. Note,
% two-tailed thresholds are used.
% Are there any permutaion resamples to perform?
Nperm = 0;

% The job split variable is how many jobs this analysis is split into for
% sending to a comuputer cluster environment. Note that the more splits
% does not mean faster computing because there is a limitationon the number
% of available cluster nodes. When the number of job splits exceeds the
% number of available jobs then jobs are placed in the queue and need to
% wait.
% If you are running this on a single (non-cluster environment) computer
% then set this to equal 1.
NJobSplit = 0;

% Specify the thresholds used in the analysis. 
% Some of the images and results in the mediation analysis perform
% parametric tests where the threshold can be more dynamically changed.
% However, for  estimation of the significance of the paths resampling
% methods are required and the probabilities need to be calculated at the
% time of the bootstrapping and therefore, made a prior.
Thresh = [0.2 0.1];

% Create a structure which will contain all information for this analysis.
ModelInfo = {};
ModelInfo.BaseDir = BaseDir;
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
ModelInfo.Nsub = size(data{1},1);
ModelInfo.Nvar = Nvar;
ModelInfo.Nvoxels = Nvoxels;

% Prepare the output data header
DataHeader.fname = '';
DataHeader.descrip = '';
DataHeader.dt = [16 0];
ModelInfo.DataHeader = DataHeader;

%% Model Specific Setup Number One
% The first model is to perform a simple mediation analysis with the
% voxel-wise brain data as the mediation between a dichotomous variable
% (Age group) and a continuous variable (Behavior). 

%     M
%    / \
%   /   \
%  X     Y

Model1 = ModelInfo;

% Name of the output folder
Tag = 'ExampleModel1_boot';
Model1.Tag = Tag;

% The first model is a basic mediation model testing whether the effect of age
% group on behavior is mediated by the voxel-wise brain measures.


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

Direct = zeros(Model1.Nvar);
Direct(1,[2 3]) = 1;
Direct(2,3) = 1;

% To specify interactions in the model create another matrix the same size
% as the DIRECT matrix.
% An interaction between two variables in predicting a third, the rows
% of the interacting variables will both have a value of one within the
% column of the variable they are predicting. Be sure to include in the
% DIRECT matrix the main effects for each of these variables on the
% predicted variable.

% For this model there are no interactions so this is kept as a zero
% matrix.
Inter = zeros(Model1.Nvar);

% The paths to be tested are included in another matrix the same size of
% the DIRECT matrix. One difference here is that the PATHS matrix may have
% a third dimension to test multiple paths within a single model. The steps
% along a path are specified with integers, i.e. step 1, step 2 ...
Paths = zeros(Model1.Nvar);
Paths(1,2) = 1;
Paths(2,3) = 2;


Model1.Direct = Direct;
Model1.Inter = Inter;
Model1.Paths = Paths;

%% Run the analysis for Model 1
% The following two functions actually run the mediation analysis and then
% writes out the results to NIFTI format images.

ResultsFolder = PrepareDataForProcess(Model1);


WriteOutResults(ResultsFolder)


%% Model Specific Setup Number 2
% The second model is to perform a simple mediation analysis with the
% voxel-wise brain data as the mediation between a dichotomous variable
% (Age group) and a continuous variable (Behavior). This analysis now
% includes covariates also.


Model2 = ModelInfo;

% Add the covariate names
Model2.Names = [Model2.Names 'Cov1' 'Cov2'];

% Add the covariate data
Model2.data = [Model2.data COVARIATES(:,1) COVARIATES(:,2)];

% Update the number of variables value 
Model2.Nvar = length(Model2.Names);

% Name of the output folder
Tag = 'ExampleModel2_boot';
Model2.Tag = Tag;

% Specify the simple mediation model as above

Direct = zeros(Model2.Nvar);
Direct(1,[2 3]) = 1;
Direct(2,3) = 1;

% Now specify the covariates in the models
Direct([4,5],2) = 1;
Direct([4,5],3) = 1;

% For this model there are no interactions so this is kept as a zero
% matrix.
Inter = zeros(Model2.Nvar);

% The paths do not change
Paths = zeros(Model2.Nvar);
Paths(1,2) = 1;
Paths(2,3) = 2;


Model2.Direct = Direct;
Model2.Inter = Inter;
Model2.Paths = Paths;
%% Run the analysis for Model 2
% The following two functions actually run the mediation analysis and then
% writes out the results to NIFTI format images.

ResultsFolder = PrepareDataForProcess(Model2);

% The writing of the images should ideally be perfomed by the cluster once
% the analyses have completed. Therefore, a check is needed or better yet a
% wait command for the cluster job: 
% e.g. qsub -hold_jid job1,job2 -N job3 ./c.sh
WriteOutResults(ResultsFolder)


%% Model Specific Setup Number Three
% The third model is to perform a moderated-mediation analysis with the
% voxel-wise brain data as the mediation between a dichotomous variable
% (Age group) and a continuous variable (Behavior). 
% Now there is an interaction between the brain measure and age group in
% predicting the behavior variable

%     M
%    / \
%   /  /\
%  /  /  \
% X --    Y

Model3 = ModelInfo;

% Update the number of variables value 
Model3.Nvar = length(Model3.Names);

% Name of the output folder
Tag = 'ExampleModel3_boot';
Model3.Tag = Tag;

% Specify the simple mediation model as above

Direct = zeros(Model3.Nvar);
Direct(1,[2 3]) = 1;
Direct(2,3) = 1;


% For this model there is interaction.
% There can actually be more then one interaction in a single model. If
% that is the case thise interaction matrix will need a third dimension. As
% of right now this is not supported.

Inter = zeros(Model3.Nvar);
Inter([1 2],3) = 1;


% The paths do not change
Paths = zeros(Model3.Nvar);
Paths(1,2) = 1;
Paths(2,3) = 2;


Model3.Direct = Direct;
Model3.Inter = Inter;
Model3.Paths = Paths;
%% Run the analysis for Model 3
% The following two functions actually run the mediation analysis and then
% writes out the results to NIFTI format images.

ResultsFolder = PrepareDataForProcess(Model3);

% The writing of the images should ideally be perfomed by the cluster once
% the analyses have completed. Therefore, a check is needed or better yet a
% wait command for the cluster job: 
% e.g. qsub -hold_jid job1,job2 -N job3 ./c.sh
WriteOutResults(ResultsFolder)

%% Perform Model 1 analyses on a single data point and print out the results 

% It is also possible to use this software to do single point analyses. The
% model setup is the same as above. The real difference is that the number
% of voxels and the list of indices need to be set to one to indicate that
% a "single voxel" is being analyzed.

SinglePointModel = Model1;
% extract a single voxel.
SinglePointModel.data{2} = Model1.data{2}(:,45);
SinglePointModel.Indices = 1;
SinglePointModel.Nvoxels = 1;
SinglePointModel.Tag = 'SinglePointModel';

% The data is entered in the above examples as a cell array of arrays or
% vectors. The first thing to do is to collapse the cell array into a
% single array for running the regression models.
SinglePointModel = ExtractDataFromVoxel(SinglePointModel,1);

% This function is at the core of the multi-voxel analyses.
Results = OneVoxelProcessBootstrap(SinglePointModel);

% This function prints the results to the screen.
PrintResults(SinglePointModel,Results)
%% Perform Model 2 analyses on a single data point and print out the results 
SinglePointModel = Model2;
SinglePointModel.data{2} = Model2.data{2}(:,45);
SinglePointModel.Indices = 1;
SinglePointModel.Nvoxels = 1;
SinglePointModel.Tag = 'SinglePointModel';
SinglePointModel = ExtractDataFromVoxel(SinglePointModel,1);

Results = OneVoxelProcessBootstrap(SinglePointModel);
PrintResults(SinglePointModel,Results)

%% Perform Model 3 analyses on a single data point and print out the results 
SinglePointModel = Model3;
SinglePointModel.data{2} = Model3.data{2}(:,45);
SinglePointModel.Indices = 1;
SinglePointModel.Nvoxels = 1;
SinglePointModel.Nboot = 10000;
SinglePointModel.Tag = 'SinglePointModel';
SinglePointModel = ExtractDataFromVoxel(SinglePointModel,1);

Results = OneVoxelProcessBootstrap(SinglePointModel);
PrintResults(SinglePointModel,Results)

%% Perform Model 1 using the cluster
% The real advantage of this software is the use of a cluster computing
% environment.
% In order to use the cluster environment you need to specify the number of
% jobs to split the analysis into.

ClusterModel1 = Model1;
% If the analysis uses bootstrapping then the data itself is broken up into
% chunks and each chunk is processed by a different cluster job. If the
% analysis uses permutation testing then the number of permutations is
% split across the different jobs.

ClusterModel1.NJobSplit = 10;
ClusterModel1.Nboot = 500;
% Using the cluster requires the use fo creating job shell scripts using
% this function:
%
%       CreateClusterJobFile(Command,fid)
%
% Once a file is created its fid (file identifier) is passed to this
% command along with the command to be executed. 
% THis function needs to be modified to make it site specific by specifying
% the install locations of MatLab and SPM.

ResultsFolder = PrepareDataForProcess(ClusterModel1);

%%
% Other models:
%     M
%     |
% X ----- Y

%    M1--M2
%   /      \
%  /        \
% X          Y

%  W  M
%   \/ \
%   /   \
%  X     Y

%     M  V
%    / \/
%   /   \
%  X     Y

%     M
%    / \
%   /\ /\
%  /  W  \
% X       Y

%    M1--M2
%   /      \
%  /       /\
% X--------  Y

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

Nvoxels = 50;
NSub = 100;
X = round(rand(NSub,1));
a = 0.5;
M = a.*repmat(X,1,Nvoxels) + randn(NSub,Nvoxels);
b = 0.5;
cP = 0.5;
Y = cP.*X + b.*M(:,4) + randn(NSub,1);

corr([X M Y]);
BaseDir = '/home/steffejr/Scripts/TestThimgs';

data = {};
data{1} = X;
data{2} = M;
data{3} = Y;

% Create a structure of corresponding names for each cell in the data
% structure
Names = {};
Names{1} = 'A';
Names{2} = 'B';
Names{3} = 'C';

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
NJobSplit = 10;

% Specify the thresholds used in the analysis. 
% Some of the images and results in the mediation analysis perform
% parametric tests where the threshold can be more dynamically changed.
% However, for  estimation of the significance of the paths resampling
% methods are required and the probabilities need to be calculated at the
% time of the bootstrapping and therefore, made a prior.
Thresh = [0.05];

% Create a structure which will contain all information for this analysis.
ModelInfo = {};
ModelInfo.BaseDir = BaseDir;
ModelInfo.Names = Names;
ModelInfo.data = data;
ModelInfo.Nboot = Nboot;
ModelInfo.Nperm = Nperm;
ModelInfo.BaseDir = BaseDir;

ModelInfo.Indices = 1;
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
ModelInfo.Nvoxels = 1;

% Prepare the output data header
DataHeader.fname = '';
DataHeader.descrip = '';
DataHeader.dt = [16 0];
ModelInfo.DataHeader = DataHeader;

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




SinglePointModel = ExtractDataFromVoxel(Model1,1);

Results = OneVoxelProcessBootstrap(SinglePointModel);
PrintResults(SinglePointModel,Results)


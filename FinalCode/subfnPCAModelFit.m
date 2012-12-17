%function subfnPCAModelFit(data)
tic
% This program needs to receive the entire data set, all voxels because it
% is doing the PCA.
%ModelNum = AnalysisParameters.ModelNum;
BaseDir = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes';
JobsDir = fullfile(BaseDir,'jobs');
DataFile = fullfile(BaseDir,'AllData')
load(DataFile)

if ~exist(JobsDir,'dir')
    % Make the outrput path for cluster jobs
    mkdir(JobsDir)
end

ModelNum = '4';
[NSub NMed NVox] = size(data.M);
if NMed > 1
    error('This doesn''t work with multiple mediators yet!');
end
NCOV = size(data.COV,2);

NPCs = 12;
% Create matrix of all possible cominations of PCs
combo_matrix = boolean_enumeration_f(NPCs);
% how many combos are there?
NCombos = size(combo_matrix,1);

remove_row_means = 1;




%% Leave one out model selection
fprintf(1,'** Starting the leave one out process **\n');

% This is the SLOWEST method that performs the PCA for every LOO and fits
% the regression model for every combination. Therefore, NO shortcuts are
% taken.
tic
subfnClusterLOOCV(DataFile, NSub, JobsDir, NPCs)
NDone = 0;
h = waitbar(NDone/NSub,'Running LOOCV ...');
while NDone < NSub
    NDone = length(dir('LOOCV*.mat'));
    waitbar(NDone/NSub,h);
    pause(1);
end
close(h)
% Find the resultant LOOCV files and load them
selected_PCs_LOOCV = subfnCompileClusterLOOCV(JobsDir,NPCs,NSub)
toc

 
%% AIC Model selection
% do the AIC model selection approach
combo_matrix = boolean_enumeration_f(NPCs);
NCombos = size(combo_matrix,1);
AICmatrix = zeros(NCombos,1);
LOOCVmatrix = zeros(NCombos,1);
[lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(data.M)', remove_row_means);
ssf = squeeze(data.M) * eigenimages_noZeroes;
ssfSubset = ssf(:,1:NPCs);

for j = 1:NCombos
    selected_PCs = find(combo_matrix(j,:));

    TrainData = ssfSubset(:,selected_PCs);
    %behav_fit_coef = FullModelbehav_fit_coef([1 selected_PCs+1 NPCs+2:end]);
    S = subfnregstats(data.Y,[TrainData data.X]);
    AICmatrix(j) = S.AIC;
    LOOCVmatrix(j) = S.CV;
end
selected_PCs_AIC = find(combo_matrix(find(AICmatrix == min(AICmatrix)),:));

%% Create the point estimate image
% perform the PCA on the original full data set

selected_PCs = [1     2     7     8     9    10    11    12]
selected_PCs = selected_PCs_LOOCV;
remove_row_means = 1;
[lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(data.M)', remove_row_means);
ssf = squeeze(data.M) * eigenimages_noZeroes;
ssfSubset = ssf(:,1:NPCs);

PE_behav_fit_coef = subfnregress(data.Y,[ssfSubset(:,selected_PCs) data.X]);
% create the SSF image
temp = eigenimages_noZeroes(:, selected_PCs) * PE_behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent

PE_behav_fit_composite_PC_image = temp / norm(temp);
clear temp;
%%%%% Obtain SSFs associated with the normalized best linear behavioral-fit image %%%%%
PE_behav_fit_composite_PC_image_ssfs = squeeze(data.M) * PE_behav_fit_composite_PC_image;

%%
tempdata = data;
tempdata.names.M={'SSFs'};
tempdata.names.Y='medRT';
tempdata.M = PE_behav_fit_composite_PC_image_ssfs;
tempdata.Thresholds = [0.05];
Parameters = subfnVoxelWiseProcessBatch(tempdata);
subfnPrintResults(Parameters{1})


[error img] = Commonality_2Pred(tempdata.X,tempdata.M,tempdata.Y,{tempdata.names.X,tempdata.names.M{1},tempdata.names.Y});

%% Now use this best set of PCs for bootstrapping
% But this needs to create the mediating pattern. Therefore, the SSF image
% needs to be weighted by the other effects in the model

Nboot = 50;
%selected_PCs = [1 2 4 5];
NTimes = 40;

subfnClusterBootStrapPC(data,Nboot,selected_PCs,NTimes,JobsDir)
NDone = 0;
h = waitbar(NDone/NTimes,'Running Bootstrapping ...');
while NDone < NTimes
    NDone = length(dir(fullfile(JobsDir,'BootStrapChunk*.mat')));
    waitbar(NDone/NTimes,h);
    pause(1);
end
close(h)
BootData = subfnCompileClusterBootStrap(JobsDir);


%%
clear eigenimages_noZeroes tempdata
% create variance map
Vmap = var(BootData,0,2);
% create standard deviation map
SDmap = std(BootData,0,2);
% Create Z-map
Zmap = PE_behav_fit_composite_PC_image./SDmap;
% create BCaci confidence intervals for each voxel

% create percentile based confidecne intervals for each voxel
toc

cd DualmSing_XageGr_MfMRI_YmedRT_10-Dec-2012_14-50
load AnalysisParameters
cd ..

Vo=AnalysisParameters.V;
Vo.n = [1 1];
Vo.descrip='';
Vo.fname=fullfile(pwd,'med_PCA_12451516_Zmap.nii');
Y = zeros(Vo.dim);
Y(AnalysisParameters.Indices)=Zmap;
spm_write_vol(Vo,Y)

% for i = 1:Nsub
%   Leave one out matrix=(Nsub-x-Ncombinations_of_PCs)
%   Apply PCA
%   for j = 1:all combinations of PCs
%       fit the regression (linear/iteratively)
%       predict left out subject
%       find squared error (one number)
%   end
% end
% sum over all subjects
% find the combination of PCs that has the lowest error
%
% for i = 1:Nboot
%   resample
%   Apply PCA
%   pick combo of PCs as calculated in the leave one out procedure
%   create the combo Image
%   save this image to the bootstrap (HUGE) resample array
%           matrix=(Nboot-x-Nvoxels)
%   fit the regression (linear/iteratively)
%   for k = 1:Nprobe
%       adjust the moderator term
%       fit the regression (linear/iteratively)
%       save parameter estimates
%   end
%   save parameter estimates
% end
%
% calculate the Z-map combo image using the bootstrap standard error
% calculate the confidence intervals for each voxel of the combo image
% calculate the confidence intervals for each regression parameter and the
% conditional indirect effect
% save to file:
%   parameter maps


%   conditional indirect effect maps at the probe values
%   confidence interval maps for the parameter maps
%   p-value maps for the parameter maps
%   Z-map for the combo image
%   combo image confidence intervals
%   combo image p-value maps




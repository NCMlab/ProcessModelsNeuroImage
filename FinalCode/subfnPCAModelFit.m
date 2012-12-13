%function subfnPCAModelFit(data)
tic
% This program needs to receive the entire data set, all voxels because it
% is doing the PCA.
%ModelNum = AnalysisParameters.ModelNum;
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

% Create bootstrap resamples
Nboot = 1000;
BootStrapResamples = zeros(NSub,Nboot,'uint16');
if ~isempty(data.STRAT)
    Gr1 = find(data.STRAT == 0);
    NGr1 = length(Gr1);
    Gr2 = find(data.STRAT == 1);
    NGr2 = length(Gr2);
else
    NGr1 = [];
    NGr2 = [];
end
for i = 1:Nboot
    if isempty(data.STRAT)
        Samp =  floor(N*rand(N,1))+1;
    else
        Samp1 = floor(NGr1*rand(NGr1,1))+1;
        Samp2 = floor(NGr2*rand(NGr2,1))+1+NGr1;
        Samp = [Samp1; Samp2];
    end
    BootStrapResamples(:,i) = Samp;
end
% Create the array that will contain all the bootstrap resample IMAGES
BootStrapResampleImages = zeros(NVox,Nboot,'single');


%% Leave one out model selection
LOOerrorMatrix = zeros(NSub, 2^NPCs-1);
tLOOerrorMatrix = zeros(NSub, 2^NPCs-1);
fprintf(1,'** Starting the leave one out process **\n');

% This is the SLOWEST method that performs the PCA for every LOO and fits
% the regression model for every combination. Therefore, NO shortcuts are
% taken.
for i = 1:NSub
    % everything in this loop should be made into a job for cluster
    % submission and the only thing returned is the LOO values. Each call
    % creates its own combo_matrix
    %
    % Inputs:
    %   data
    %   left out subject
    %   NPCs
    % Outputs: LOO vector
    fprintf(1,' working on subject %4d of %4d, ',i,NSub);
    % Create leave one out matrix
    CurrentSubjects = [1:NSub];
    CurrentSubjects(i) = 0;
    CurrentSubjects = find(CurrentSubjects);
    TrainData = [];
    TrainData.X = data.X(CurrentSubjects);
    TrainData.M = squeeze(data.M(CurrentSubjects,:,:));
    TrainData.Y = data.Y(CurrentSubjects);
    if NCOV
        TrainData.COV = data.COV(CurrentSubjects);
    else
        TrainData.COV = [];
    end
    TrainData.ModelNum = ModelNum;
    % apply PCA
    % but squeeze out the multiple mediator dimension
    [lambdas, eigenimages_noZeroes, w] = pca_f(TrainData.M', remove_row_means);
    TrainSSF = squeeze(TrainData.M) * eigenimages_noZeroes;
    TrainSSFsubset = TrainSSF(:,1:NPCs);
    % cycle over all combinations of PCs
    testdata = TrainData;
    testdata.M = TrainSSFsubset;
    b = subfnCallRegressPCs(testdata);
    tic
    for j = 1:NCombos
        selected_PCs = find(combo_matrix(j,:));
        NPCsSubset = length(selected_PCs);
        b_subset = [b(1); b(1+selected_PCs); b(NPCs+2:end)];
        tLOOerrorMatrix(i,j) = (data.Y(i)...
            - squeeze(data.M(i,:,:))'*eigenimages_noZeroes(:,selected_PCs)*b_subset(2:1+length(selected_PCs))...
            - data.X(i)*b_subset(1))^2;
    end
    t1 = toc;
    fprintf(1,'%0.3f sec, ',t1)
    tic
    for j = 1:NCombos
        selected_PCs = find(combo_matrix(j,:));
        % create a copy of the data with the left out subject so that the
        % voxel-wise data can be replaced with the SSFs of the PCs
        LOOdata = TrainData;
        LOOdata.M = TrainSSFsubset(:,selected_PCs);
        behav_fit_coef = subfnCallRegressPCs(LOOdata);
        LOOerrorMatrix(i,j) = (data.Y(i)...
            - squeeze(data.M(i,:,:))'*eigenimages_noZeroes(:,selected_PCs)*behav_fit_coef(2:1+length(selected_PCs))...
            - data.X(i)*behav_fit_coef(1))^2;
    end
     t2 = toc;
    fprintf(1,'%0.3f sec\n',t2)
end

%%
%sLOOerrorMatrix1 = sum(LOOerrorMatrix);
sLOOerrorMatrix = sum(LOOerrorMatrix);
stLOOerrorMatrix = sum(tLOOerrorMatrix);
selected_PCs_LOO = find(combo_matrix(find(sLOOerrorMatrix == min(sLOOerrorMatrix)),:));
tselected_PCs_LOO = find(combo_matrix(find(stLOOerrorMatrix == min(stLOOerrorMatrix)),:));


% AIC Model selection
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
selected_PCs_LOOCV = find(combo_matrix(find(LOOCVmatrix == min(LOOCVmatrix)),:));
[sLOOerrorMatrix' AICmatrix LOOCVmatrix]


%% Create the point estimate image
% perform the PCA on the original full data set

selected_PCs = selected_PCs_LOO
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

%% Now use this best set of PCs for bootstrapping
% But this needs to create the mediating pattern. Therefore, the SSF image
% needs to be weighted by the other effects in the model
for i = 1:200
    clear tempdata
    TrainData.X = data.X(BootStrapResamples(:,i));
    TrainData.M = data.M(BootStrapResamples(:,i),:,:);
    TrainData.Y = data.Y(BootStrapResamples(:,i));
    TrainData.ModelNum = ModelNum;
    % Apply PCA
    [lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(TrainData.M)', remove_row_means);
    % calculate the subject scaling factors
    ssf = squeeze(TrainData.M) * eigenimages_noZeroes;
    ssfSubset = ssf(:,1:NPCs);
    behav_fit_coef = subfnregress(TrainData.Y, [ssfSubset(:,selected_PCs) TrainData.X]);
    % create the SSF image
    temp = eigenimages_noZeroes(:, selected_PCs) * behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent
    % Created the weighted sum pf SSFs regressor 
    M = ssfSubset(:,selected_PCs)*behav_fit_coef(1 + 1:1 + length(selected_PCs));
    a_coef = subfnregress(M,TrainData.X);
    % Calculate the mediated effect
    temp = temp.*a_coef(end);
    behav_fit_composite_PC_image = temp / norm(temp);
    clear temp;
    BootStrapResampleImages(:,i) = behav_fit_composite_PC_image;
    if ~mod(i,100)
        fprintf(1,'Bootstrap resample: %5d\n',i);
    end
end
%%
clear eigenimages_noZeroes tempdata
% create variance map
Vmap = var(BootStrapResampleImages,0,2);
% create standard deviation map
SDmap = std(BootStrapResampleImages,0,2);
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
Vo.fname=fullfile(pwd,'med_PCA_1245_Zmap.nii');
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




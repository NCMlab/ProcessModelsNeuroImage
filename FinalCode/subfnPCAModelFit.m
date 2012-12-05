%function subfnPCAModelFit(data)
tic
% This program needs to receive the entire data set, all voxels because it
% is doing the PCA.
ModelNum = AnalysisParameters.ModelNum;
ModelNum = '4';
[NSub NMed NVox] = size(data.M);
if NMed > 1
    error('This doesn''t work with multiple mediators yet!');
end
NPCs = 6;
% Create matrix of all possible cominations of PCs
combo_matrix = boolean_enumeration_f(NPCs);
% how many combos are there?
NCombos = size(combo_matrix,1);
remove_row_means = 1;
LOOerrorMatrix = zeros(NSub,NCombos);
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
fprintf(1,'** Starting the leave one out process **\n');    
for i = 1:NSub
    tic
    % everything in this loop should be made into a job for cluster
    % submission and the only thing returned is the LOO values. Each call
    % creates its own combo_matrix
    % 
    % Inputs: 
    %   data
    %   left out subject
    %   NPCs
    % Outputs: LOO vector
    fprintf(1,' working on subject %4d of %4d: ',i,NSub);
    % Create leave one out matrix
    CurrentSubjects = [1:NSub];
    CurrentSubjects(i) = 0;
    CurrentSubjects = find(CurrentSubjects);
    tempdata.X = data.X(CurrentSubjects);
    tempdata.M = data.M(CurrentSubjects,:,:);
    tempdata.Y = data.Y(CurrentSubjects); 
    tempdata.ModelNum = ModelNum;
    % apply PCA
    % but squeeze out the multiple mediator dimension
    [lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(tempdata.M)', remove_row_means);
    ssf = squeeze(tempdata.M) * eigenimages_noZeroes;
    ssfSubset = ssf(:,1:NPCs);
    
    for j = 1:NCombos
            % Select the current combination os PCs
        if AnalysisParameters.ModelNum == '4'
            % This is the easy case where simple linear regression can be
            % used instead of an iterative model fit.
            % Used the selected combination of SSFs for this
            selected_PCs = find(combo_matrix(j,:));
            behav_fit_coef = subfnregress(tempdata.Y, [ssfSubset(:,selected_PCs) tempdata.X]);
            % create the SSF image
            temp = eigenimages_noZeroes(:, selected_PCs) * behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent
            behav_fit_composite_PC_image = zeros(NVox, 1);
            behav_fit_composite_PC_image = temp / norm(temp);
            %%%%% Obtain SSFs associated with the normalized best linear behavioral-fit image %%%%%
            behav_fit_composite_PC_image_ssfs = squeeze(tempdata.M) * behav_fit_composite_PC_image;
            % forward apply this SSF to the left out subjects raw data
            % predict the left out subject
            predictedValue = squeeze(data.M(i,:,:))'*behav_fit_composite_PC_image;
            LOOerrorMatrix(i,j) = [predictedValue + behav_fit_coef(1) - data.Y(i)]^2;
        end
    end
    t = toc;
    fprintf(1,'%0.2f sec\n',t);
end
sLOOerrorMatrix = sum(LOOerrorMatrix);
selected_PCs = find(combo_matrix(find(sLOOerrorMatrix == min(sLOOerrorMatrix)),:));
fprintf(1,'The optimal PCs were selected\n');
%% Create the point estimate image
% perform the PCA on the original full data set
[lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(data.M)', remove_row_means);
ssf = squeeze(data.M) * eigenimages_noZeroes;
ssfSubset = ssf(:,1:NPCs);

PE_behav_fit_coef = subfnregress(data.Y,[ssfSubset(:,selected_PCs) data.X]);
% create the SSF image
temp = eigenimages_noZeroes(:, selected_PCs) * PE_behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent
PE_behav_fit_composite_PC_image = zeros(NVox, 1);
PE_behav_fit_composite_PC_image = temp / norm(temp);
clear temp;
%%%%% Obtain SSFs associated with the normalized best linear behavioral-fit image %%%%%
PE_behav_fit_composite_PC_image_ssfs = squeeze(data.M) * PE_behav_fit_composite_PC_image;

%% Now use this best set of PCs for bootstrapping
for i = 1:Nboot
    clear tempdata
    tempdata.X = data.X(BootStrapResamples(:,i));
    tempdata.M = data.M(BootStrapResamples(:,i),:,:);
    tempdata.Y = data.Y(BootStrapResamples(:,i));
    tempdata.ModelNum = data.ModelNum;
    % Apply PCA
    [lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(tempdata.M)', remove_row_means);
    % calculate the subject scaling factors
    ssf = squeeze(tempdata.M) * eigenimages_noZeroes;
    ssfSubset = ssf(:,1:NPCs);
    behav_fit_coef = subfnregress(tempdata.Y, [ssfSubset(:,selected_PCs) tempdata.X]);
    % create the SSF image
    temp = eigenimages_noZeroes(:, selected_PCs) * behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent
    behav_fit_composite_PC_image = zeros(NVox, 1);
    behav_fit_composite_PC_image = temp / norm(temp);
    clear temp;
    BootStrapResampleImages(:,i) = behav_fit_composite_PC_image;
end
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




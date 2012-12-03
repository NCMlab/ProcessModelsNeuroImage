function subfnPCAModelFit(data)
% This program needs to receive the entire data set, all voxels because it
% is doing the PCA.

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
for i = 1:NSub
    % Create leave one out matrix
    CurrentSubjects = [1:NSub];
    CurrentSubjects(i) = 0;
    CurrentSubjects = find(CurrentSubjects);
    tempdata.X = data.X(CurrentSubjects);
    tempdata.M = data.M(CurrentSubjects,:,:);
    tempdata.Y = data.Y(CurrentSubjects); 
    tempdata.ModelNum = data.ModelNum;
    % apply PCA
    % but squeeze out the multiple mediator dimension
    [lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(tempdata.M)', remove_row_means);
    ssf = squeeze(tempdata.M) * eigenimages_noZeroes;
    ssfSubset = ssf(:,1:NPCs);
    for j = 1:NCombos
            % Select the current combination os PCs
        if tempdata.ModelNum == '4'
            % This is the easy case where simple linear regression can be
            % used instead of an iterative model fit.
            % Used the selected combination of SSFs for this
            selected_PCs = find(combo_matrix(j,:));
            behav_fit_coef = subfnregress(tempdata.Y,[ssfSubset(:,selected_PCs)]);
            % create the SSF image
            temp = eigenimages_noZeroes(:, selected_PCs) * behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent
            behav_fit_composite_PC_image              = zeros(NVox, 1);
            behav_fit_composite_PC_image = temp / norm(temp);
            clear temp;
            %%%%% Obtain SSFs associated with the normalized best linear behavioral-fit image %%%%%
            behav_fit_composite_PC_image_ssfs = squeeze(tempdata.M) * behav_fit_composite_PC_image;
            % forward apply this SSF to the left out subjects raw data
            % predict the left out subject
            predictedValue = squeeze(data.M(i,:,:))'*behav_fit_composite_PC_image;
            LOOerrorMatrix(i,j) = [predictedValue + behav_fit_coef(1) - data.Y(i)]^2;
        end
    end
end
sLOOerrorMatrix = sum(LOOerrorMatrix);
best_set = combo_matrix(find(sLOOerrorMatrix == min(sLOOerrorMatrix)),:);

% Now use this best set of PCs 
for i = 1:Nboot
end

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




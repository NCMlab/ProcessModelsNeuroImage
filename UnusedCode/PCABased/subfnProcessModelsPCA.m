function subfnProcessModelsPCA(data)
% data needs to contain the original data
[NSub NPC] = size(data.M);
% Create the set of all possible PC combinations

power_set = boolean_enumeration_f(NPC);
NCombos = size(power_set,1)
LOOCVmatrix = zeros(NSub,NCombs);
remove_row_means = 1;
for i = 1:Nsub
%   Leave one out matrix=(Nsub-x-Ncombinations_of_PCs)
for j = 1:NCombos
    tempSubjectList = 1:NSub;
    tempSubjectList(i) = 0;
    tempSubjectList = find(tempSubjectList)
    target_images = data.M(tempSubjectList,1,:);
    meanful_set = find(power_set(j,:));
    [lambdas, eigenimages_complete, w] = pca_f(target_images, remove_row_means, meanful_set)
%   Apply PCA
%   for j = 1:all combinations of PCs
%       fit the regression (linear/iteratively)
%       predict left out subject
%       find squared error (one number)
%   end
% end
% sum over all subjects
% find the comboination of PCs that has the lowest error
%
% for i = 1:Nboot
%   resample
%   Apply PCA
%   pick combo of PCs as calculated in the leave one out procedure
%   create the combo Image
%   save this image to the bootstrap (HUGE) resample array
%           matrix=(Nboot-x-Nvoxels)
%   fit the regression (linear/iteratively)
%       iterative function: F = subfnRegressPCs(coef,data)
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
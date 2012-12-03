function [lambdas, eigenimages_complete, w] = pca_f(target_images, remove_row_means, meanful_set)
%
%function [lambdas, eigenimages_complete, w] = pca_f(target_images, remove_row_means, meanful_set)
%
%%%%% Performs principal component analysis (without any design matrix).
%%%%% Target image array is doubly-normalized (with both row and column means removed).
%%%%% "target_images"       : target image array, with voxels across subjects/tasks in rows
%%%%% "meanful_set"         : location indices of voxels(rows) with zero entries
%%%%% "remove_row_means"    : remove mean image (row means in the target image array) from all
%%%%%                         images ("0": no removal; "1": removal)
%%%%% "lambdas"             : eigenvalues in descending order, without normalization to unity
%%%%% "eigenimages_complete": normalized eigenimages, represented in columns, with the same
%%%%%                         number of voxels(rows) as that of the target image array
%%%%% "w"                   : eigenvectors that correspond to eigenvalues("lambdas")
%%%%% Note: As a consequence of removing row means, the number of eigenvalues/eigenimages is 
%%%%%       one fewer than the number of target images.



if (remove_row_means ~= 0) & (remove_row_means ~= 1)
   warning('The input argument "remove_row_means" has to hold either the value "0" or the value "1".')
   disp(' ')
   lambdas = [];
   eigenimages_complete = [];
   w = [];
   return;
end

if nargin == 2
   temp = sum(target_images, 2);
   meanful_set = find(temp ~= 0);
   clear temp;
end



data = target_images(meanful_set, :);	%%%%% Remove voxels with zero activation values %%%%%

data_matrix_dimensions = size(data);

if remove_row_means == 1
   %%%%% Remove mean brain image activation for each voxel and for each subject/task %%%%%
   %%%%% Each element of any given row represents the same voxel across subject/task %%%%%
   %%%%% Each element of any given column represents a distinct voxel in a single    %%%%%
   %%%%%   subject/task                                                              %%%%%
   %%%%% Note: Measuring task and subject effects in voxel covariance pattern.       %%%%%
   mddata = remove_row_column_means_f(data);
else
   %%%%% Note: Measuring only mean effect in voxel covariance pattern. %%%%%
   mddata = remove_column_mean_f(data);
end



%%%%% Note: With only column means removed, "cov_mat" should equal "cov(data)". %%%%%
cov_mat = mddata' * mddata / (data_matrix_dimensions(1) - 1);

%%%%% Perform eigen decomposition                                                            %%%%%
%%%%% Outputs("w") are sorted according to each component's variance contribution("lambdas") %%%%%
[w, lambdas] = eig_habeck(cov_mat);

%%%%% Obtain normalized eigenimages; each column represents one eigenimage; exclude %%%%%
%%%%%   zero eigenimages                                                              %%%%%
temp = mddata * w;
nOfPCs = rank(temp);
for i=1:nOfPCs
   temp(:, i) = temp(:, i) / norm(temp(:, i));
end
eigenimages_complete                 = zeros(size(target_images, 1), nOfPCs);
eigenimages_complete(meanful_set, :) = temp(:, 1:nOfPCs);

%%%%% Remove zero eigenvalues %%%%%
lambdas = lambdas(1:nOfPCs);

%%%%% Output warning if there is an unusual number of eigenimages %%%%%
if ( (remove_row_means == 1) & (nOfPCs ~= size(target_images, 2) - 1) ) ...
   | ...
   ( (remove_row_means == 0) & (nOfPCs ~= size(target_images, 2)) )

   disp(['     Note: The number of non-zero eigenvalues is ', num2str(nOfPCs), '.'])
end

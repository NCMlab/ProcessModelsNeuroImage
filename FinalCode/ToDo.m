% To do
% 
% Some way to calculate the clusters using SPM.
%
% Additions need to be to allow for X or Y to be voxelwise data instead of
% only M.
% The program: subfnVoxelWiseProcessBatch needs to be split into 3 parts.
% 1) make sure all the parameters are set then cycles over voxels and calls the data preparation
% and stats  = subfnVoxelWiseProcessBatch
% This program also needs to prepare the data to avoid passing large
% amounts of data between programs.
% 2) actually does all the stats = subfnVoxelWiseProcess
%
% Need to allow for voxelwise: X, M, Y, V, W, covariates
% Add the number of mediators, subjects and voxels to the data structure so
% the program does not need to figure these out. 
% For voxel the program needs to pull out the data for that voxel if that
% variable is voxelwise data. This can be a series of if statements for
% each variable for each voxel.
% In the setup set flags for each variable based on whether they are
% voxelwise or not. Then when the data fort a single voxel is being
% prepared there is just a boolean check is Xvoxelwise, then data.X +
% data(:,:,i); else X=X.

% 
% subfnVoxelWiseProcessBatch
%  This program decides whether the input is a structure or a file path.
%  If it is a structure then the data from this single point is processed.
%  If it is a file path then it is assumed that this is a "chunk" of data from a full image.
%  This program cycles over voxels, if needed, and formats the data for statistical analysis 
%      by the program:
%      subfnProcess
%      
%      
% subfnRunVoxelwiseProcess
%  This program takes analysis parameters and data and prepares it to be analyzed 
%  by the cluster computer. Therefore, it splits the data into chunks and saves the data chunks. 
%  Then each chunk is sent to the program: subfnVoxelWiseProcessBatch and processed by 
%  the cluster.
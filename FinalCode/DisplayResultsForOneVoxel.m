function [Parameters1 data1 SelectedPath] = DisplayResultsForOneVoxel(VOXmm,SelectedPath)
if nargin < 2
    SelectedPath = spm_select(1,'dir');

end
cd(SelectedPath)
if exist('AnalysisParameters.mat')
    load AnalysisParameters
else
    errordlg('this folder does not have the required AnalysticParameters.mat file');
end
% Get the header image from one image
F = dir(fullfile(SelectedPath,'*.nii'));


V = AnalysisParameters.V;
%VOXmm = [30 56 8];
% Convert to locations
VOXind = inv(V.mat)*[VOXmm 1]';
% Find the index
[index] = sub2ind(V.dim,VOXind(1),VOXind(2),VOXind(3));
% Find the data for this index
NVox = AnalysisParameters.Nvoxels;
NVoxPerChunk = ceil(NVox/AnalysisParameters.NJobSplit);
DataIndex = find(AnalysisParameters.Indices == index);
DataFileIndex = ceil(DataIndex/NVoxPerChunk);
% Load the data file
load(fullfile(SelectedPath,[sprintf('data_%04d',DataFileIndex)]))
% Load the results file
load(fullfile(SelectedPath,[sprintf('Results_%04d',DataFileIndex)]))
% Find the indeex in the data file
dataChunkIndex = find(data.Indices == index);
data1 = data;
Parameters1 = Parameters{dataChunkIndex};
data1.M = data.M(:,:,dataChunkIndex);
% Reprocess the data location
%Parameters = subfnVoxelWiseProcessBatch(data1);
% Print the results to the screen

subfnPrintResults(Parameters1)
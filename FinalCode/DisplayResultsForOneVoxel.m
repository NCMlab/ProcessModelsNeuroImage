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

%%
groups = unique(data1.X);
if length(groups) == 2
    % dicotomous, print mean differences
    fprintf(1,'\n**** GROUP DIFFERENCES IN M ****\n');
    Gr1 = find(data1.X == groups(1));
    Gr2 = find(data1.X == groups(2));
    [Mh Mp Mci Mstats] = ttest2(data1.M(Gr1),data1.M(Gr2));
    [Mh1 Mp1 Mci1 Mstats1] = ttest(data1.M(Gr1));
    [Mh2 Mp2 Mci2 Mstats2] = ttest(data1.M(Gr2));
    
    M1 = mean(data1.M(Gr1));
    M2 = mean(data1.M(Gr2));
    S1 = std(data1.M(Gr1));
    S2 = std(data1.M(Gr2));
    % print this all
    fprintf(1,'%10s%10s%10s%10s%10s\n','group','mean','std','t-val','p-val');
    fprintf(1,'%10s%10.4f%10.4f%10.4f%10.4f\n','Gr1',M1,S1,Mstats1.tstat,Mp1);
    fprintf(1,'%10s%10.4f%10.4f%10.4f%10.4f\n','Gr2',M2,S2,Mstats2.tstat,Mp2);    
    fprintf(1,'%10s%10s%10.4f,%10s%10.4f\n','Between:','t-val',Mstats.tstat,'p-val',Mp);
    
    [rA pA] = corr([data1.M data1.Y]);
    
    [r1 p1] = corr([data1.M(Gr1) data1.Y(Gr1)]);
    [r2 p2] = corr([data1.M(Gr2) data1.Y(Gr2)]);
    fprintf(1,'corr M-Y group1: %0.4f, p-val %0.4f\n',r1(1,2),p1(1,2));
    fprintf(1,'corr M-Y group2: %0.4f, p-val %0.4f\n',r2(1,2),p2(1,2));
    fprintf(1,'corr M-Y All: %0.4f, p-val %0.4f\n',rA(1,2),pA(1,2));
    
    
end
    


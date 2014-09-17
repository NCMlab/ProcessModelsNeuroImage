function [Parameters1 SinglePointModel SelectedPath] = DisplayResultsForOneVoxel(VOXmm,SelectedPath,PrintFlag)
if nargin < 2
    SelectedPath = spm_select(1,'dir');
    PrintFlag = 1;
elseif nargin == 2
    PrintFlag = 1;
end
cd(SelectedPath)
cd('data')
if exist('ModelInfo.mat')
    load ModelInfo
else
    errordlg('this folder does not have the required AnalysticParameters.mat file');
end


% Get the header information
V = ModelInfo.DataHeader;

% Convert input mm location to voxel location
VOXind = inv(V.mat)*[VOXmm 1]';

% Find the index
[index] = sub2ind(V.dim,VOXind(1),VOXind(2),VOXind(3));
fprintf(1,'Using index: %d\n',index); 

% Find the data for this index
% The data format will differ based on the type of analysis. 
% Therefore, determine if this is a bootstrap or a permutation test.
if ModelInfo.Nboot > 0 
    ModelType = 'bootstrap';
elseif ModelInfo.Nperm > 0 
    ModelType = 'permutation';
else
    error('Unknown results!');
end

NVox = ModelInfo.Nvoxels;

DataIndex = find(ModelInfo.Indices == index);
if isempty(DataIndex)
    DataIndex = 1;
end

% Is the data broken into chunks or not
if ModelInfo.NJobSplit > 1
    DataSplitFlag = 1;
else
    DataSplitFlag = 0;
end
switch ModelType
    case 'bootstrap'
        if DataSplitFlag
            NVoxPerChunk = ceil(NVox/ModelInfo.NJobSplit);
            DataFileIndex = ceil(DataIndex/NVoxPerChunk);
            
            % Load the data file
            load(fullfile(SelectedPath,'data',[sprintf('data_%04d',DataFileIndex)]))
            
            % Load the results file
            load(fullfile(SelectedPath,'Results',[sprintf('BootStrapResults_%04d',DataFileIndex)]))
            
            
            % Find the index in the data file
            dataChunkIndex = find(subModelInfo.Indices == DataIndex);
            data1 = subModelInfo;
            
            
            Parameters1 = Results{dataChunkIndex};
            
            SinglePointModel = ExtractDataFromVoxel(subModelInfo,dataChunkIndex);
        else
            DataFileIndex = 1;
            
            % Load the results file
            F = dir(fullfile(SelectedPath,'Results',sprintf('BootStrap*')));
            load(fullfile(SelectedPath,'Results',F.name));
            Parameters1 = Parameters{DataIndex};
            SinglePointModel = ExtractDataFromVoxel(ModelInfo,DataIndex);
        end
end
% Reprocess the data location
%Parameters = subfnVoxelWiseProcessBatch(data1);
% Print the results to the screen
if PrintFlag == 1
    PrintResults(SinglePointModel,Parameters1)
end



%%
% if PrintFlag == 1
%     groups = unique(data1.X);
%     if length(groups) == 2
%         % dicotomous, print mean differences
%         fprintf(1,'\n**** GROUP DIFFERENCES IN M ****\n');
%         Gr1 = find(data1.X == groups(1));
%         Gr2 = find(data1.X == groups(2));
%         [Mh Mp Mci Mstats] = ttest2(data1.M(Gr1),data1.M(Gr2));
%         [Mh1 Mp1 Mci1 Mstats1] = ttest(data1.M(Gr1));
%         [Mh2 Mp2 Mci2 Mstats2] = ttest(data1.M(Gr2));
%         
%         M1 = mean(data1.M(Gr1));
%         M2 = mean(data1.M(Gr2));
%         S1 = std(data1.M(Gr1));
%         S2 = std(data1.M(Gr2));
%         % print this all
%         fprintf(1,'%10s%10s%10s%10s%10s\n','group','mean','std','t-val','p-val');
%         fprintf(1,'%10s%10.4f%10.4f%10.4f%10.4f\n','Gr1',M1,S1,Mstats1.tstat,Mp1);
%         fprintf(1,'%10s%10.4f%10.4f%10.4f%10.4f\n','Gr2',M2,S2,Mstats2.tstat,Mp2);
%         fprintf(1,'%10s%10s%10.4f,%10s%10.4f\n','Between:','t-val',Mstats.tstat,'p-val',Mp);
%         
%         [rA pA] = corr([data1.M data1.Y]);
%         
%         [r1 p1] = corr([data1.M(Gr1) data1.Y(Gr1)]);
%         [r2 p2] = corr([data1.M(Gr2) data1.Y(Gr2)]);
%         fprintf(1,'corr M-Y group1: %0.4f, p-val %0.4f\n',r1(1,2),p1(1,2));
%         fprintf(1,'corr M-Y group2: %0.4f, p-val %0.4f\n',r2(1,2),p2(1,2));
%         fprintf(1,'corr M-Y All: %0.4f, p-val %0.4f\n',rA(1,2),pA(1,2));
%         
%         
%     end
%     
% end

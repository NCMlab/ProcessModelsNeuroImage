function subfnRunVoxelwiseProcess(AllData,AnalysisParameters)

%addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode

%addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode
%addpath /Users/jason/Desktop/ProcessModelsNeuroImage/FinalCode

BaseDir = AnalysisParameters.BaseDir;
Nsub = AnalysisParameters.Nsub;
Nmed = AnalysisParameters.Nmed;
Nvoxels = AnalysisParameters.Nvoxels;
NJobSplit = AnalysisParameters.NJobSplit;
ModelNum = AllData.ModelNum;


% how many jobs to split this into
NvoxelsPerJob = ceil(Nvoxels/NJobSplit);
% create a folder for data and results
c = clock;
OutFolderName = sprintf('%s_%s_%02d-%02d',AnalysisParameters.Tag,date,c(4),c(5))
OutFolder = fullfile(BaseDir,OutFolderName);
JobOutputFolder = fullfile(OutFolder,'JobOutput');
InJobFolder = fullfile(OutFolder,'JobFiles');
if ~exist(OutFolder)
    mkdir(OutFolder);
end
% create a folder for all the cluster job files
if ~exist(JobOutputFolder)
    mkdir(JobOutputFolder);
end
if ~exist(InJobFolder)
    mkdir(InJobFolder);
end
DataFolder = fullfile(OutFolder,'data');
if ~exist(DataFolder)
    mkdir(DataFolder);
end


% save the analysis parameters
Str = ['save ' fullfile(OutFolder,'AnalysisParameters') ' AnalysisParameters'];
eval(Str);

AllParameters = cell(Nvoxels,1);

for i = 1:NJobSplit
    % This just makes sure that the loop is entered at least once
    if (NJobSplit > 1) & (i == NJobSplit)
        break
    end
    VoxelForThisJob = [(i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob];
    data = AllData;
    data.Nboot = data.Nboot;
    
    data = CreateDataChunk(data,VoxelForThisJob,ModelNum);
    
    data.Indices = AllData.Indices(VoxelForThisJob);
    %Parameters = subfnVoxelWiseProcessBatch(temp,ModelNum,Nboot,Thresholds);
    InTag = sprintf('data_%04d',i);
    InDataPath = fullfile(DataFolder,InTag);
    Str = ['save ' InDataPath ' data  '];
    eval(Str);
    
    % If there is only ONE job specified then DO NOT send this to the
    % cluster
    if NJobSplit > 1
        jobPath = fullfile(InJobFolder,sprintf('job_%04d.sh',i));
        fid = fopen(jobPath,'w');
        fprintf(fid,'#!/bin/bash\n');
        fprintf(fid,'#PBS -N Matlab\n');
        fprintf(fid,'#PBS -m be\n');
        fprintf(fid,'#PBS -j oe\n');
        fprintf(fid,'cd %s\n',BaseDir);
        fprintf(fid,'/usr/local/matlab/bin/matlab -nodisplay << EOF\n');
        fprintf(fid,'%s\n','addpath /share/data/data5/spm8');
        fprintf(fid,'%s\n','addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode');
        fprintf(fid,'%s\n',['subfnVoxelWiseProcessBatch(''' InDataPath ''');']);
        fprintf(fid,'exit\n');
        fprintf(fid,'EOF\n');
        fclose(fid);
        %    Str = ['! qsub  ' jobPath];
        Str = ['! qsub -q veryshort.q -p -10 -e ' JobOutputFolder ' -o ' JobOutputFolder ' -l mem_free=500M ' jobPath];
      %  unix(Str);
    else
        subfnVoxelWiseProcessBatch(InDataPath);
    end
end
if NJobSplit > 1
    % now run the last chunk of data
    VoxelForThisJob =[(i-1)*NvoxelsPerJob + 1:Nvoxels];
    data = AllData;
    
    data = CreateDataChunk(AllData,VoxelForThisJob,ModelNum);
    
    % keep track of which voxels are being processed
    data.Indices = AllData.Indices(VoxelForThisJob);
    %Parameters = subfnVoxelWiseProcessBatch(temp,ModelNum,Nboot,Thresholds);
    InTag = sprintf('data_%04d',i);
    InDataPath = fullfile(OutFolder,InTag);
    Str = ['save ' InDataPath ' data '];
    eval(Str);
    
    jobPath = fullfile(InJobFolder,sprintf('job_%04d.sh',i));
    fid = fopen(jobPath,'w');
    fprintf(fid,'#!/bin/bash\n');
    fprintf(fid,'#PBS -N Matlab\n');
    fprintf(fid,'#PBS -m be\n');
    fprintf(fid,'#PBS -j oe\n');
    fprintf(fid,'cd %s\n',BaseDir);
    fprintf(fid,'/usr/local/matlab/bin/matlab -nodisplay << EOF\n');
    fprintf(fid,'%s\n','addpath /share/data/data5/spm8');
    fprintf(fid,'%s\n','addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode');
    fprintf(fid,'%s\n',['subfnVoxelWiseProcessBatch(''' InDataPath ''');']);
    fprintf(fid,'exit\n');
    fprintf(fid,'EOF\n');
    fclose(fid);
    Str = ['! qsub -q veryshort.q -e ' JobOutputFolder ' -o ' JobOutputFolder ' -l mem_free=500M ' jobPath];
   % unix(Str);
end
%
% % Once it is all done put the data back together
% AllDoneFlag = 0;
% if ~length(dir('Results_*.mat'))<NJobSplit
%     fprintf(1,'All Done!\n');
% end
% for i = 1:NJobSplit - 1
%     clear Parameters
%     load(fullfile(BaseDir,sprintf('Results_%04d',i)))
%     AllParameters((i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob) = Parameters;
% end
% clear Parameters
% load(fullfile(BaseDir,sprintf('Results_%04d',i+1)))
% AllParameters((i)*NvoxelsPerJob + 1:Nvoxels) = Parameters;

function data = CreateDataChunk(data,VoxelForThisJob,ModelNum)
switch ModelNum
    case '1'
        if size(data.X,3) > 1
            data.X = data.X(:,:,VoxelForThisJob);
        else
            data.X = data.X;
        end
        % Check the M variable
        if size(data.M,3) > 1
            data.M = data.M(:,:,VoxelForThisJob);
        else
            data.M = data.M;
        end
        % Check the Y variable
        if size(data.Y,3) > 1
            data.Y = data.Y(:,:,VoxelForThisJob);
        else
            data.Y = data.Y;
        end
        % Check the covariates
        if size(data.COV,3) > 1
            data.COV = data.COV(:,VoxelForThisJob);
        else
            data.COV = data.COV;
        end
    case '4'
        if size(data.X,3) > 1
            data.X = data.X(:,:,VoxelForThisJob);
        else
            data.X = data.X;
        end
        % Check the M variable
        if size(data.M,3) > 1
            data.M = data.M(:,:,VoxelForThisJob);
        else
            data.M = data.M;
        end
        % Check the Y variable
        if size(data.Y,3) > 1
            data.Y = data.Y(:,:,VoxelForThisJob);
        else
            data.Y = data.Y;
        end
        % Check the covariates
        if size(data.COV,3) > 1
            data.COV = data.COV(:,VoxelForThisJob);
        else
            data.COV = data.COV;
        end
    case '6'
        if size(data.X,3) > 1
            data.X = data.X(:,:,VoxelForThisJob);
        else
            data.X = data.X;
        end
        % Check the M variable
        if size(data.M,3) > 1
            data.M = data.M(:,:,VoxelForThisJob);
        else
            data.M = data.M;
        end
        % Check the Y variable
        if size(data.Y,3) > 1
            data.Y = data.Y(:,:,VoxelForThisJob);
        else
            data.Y = data.Y;
        end
        % Check the covariates
        if size(data.COV,3) > 1
            data.COV = data.COV(:,VoxelForThisJob);
        else
            data.COV = data.COV;
        end
    case '7'
        if size(data.X,3) > 1
            data.X = data.X(:,:,VoxelForThisJob);
        else
            data.X = data.X;
        end
        % Check the M variable
        if size(data.M,3) > 1
            data.M = data.M(:,:,VoxelForThisJob);
        else
            data.M = data.M;
        end
        % Check the Y variable
        if size(data.Y,3) > 1
            data.Y = data.Y(:,:,VoxelForThisJob);
        else
            data.Y = data.Y;
        end
        % Check the W variable
        if size(data.W,3) > 1
            data.W = data.W(:,:,VoxelForThisJob);
        else
            data.W = data.W;
        end
        % Check the covariates
        if size(data.COV,3) > 1
            data.COV = data.COV(:,VoxelForThisJob);
        else
            data.COV = data.COV;
        end
        
    case '14'
        if size(data.X,3) > 1
            data.X = data.X(:,:,VoxelForThisJob);
        else
            data.X = data.X;
        end
        % Check the M variable
        if size(data.M,3) > 1
            data.M = data.M(:,:,VoxelForThisJob);
        else
            data.M = data.M;
        end
        % Check the Y variable
        if size(data.Y,3) > 1
            data.Y = data.Y(:,:,VoxelForThisJob);
        else
            data.Y = data.Y;
        end
        % Check the V variable
        if size(data.V,3) > 1
            data.V = data.V(:,:,VoxelForThisJob);
        else
            data.V = data.V;
        end
        % Check the covariates
        if size(data.COV,3) > 1
            data.COV = data.COV(:,VoxelForThisJob);
        else
            data.COV = data.COV;
        end
        
    case '58'
        if size(data.X,3) > 1
            data.X = data.X(:,:,VoxelForThisJob);
        else
            data.X = data.X;
        end
        % Check the M variable
        if size(data.M,3) > 1
            data.M = data.M(:,:,VoxelForThisJob);
        else
            data.M = data.M;
        end
        % Check the Y variable
        if size(data.Y,3) > 1
            data.Y = data.Y(:,:,VoxelForThisJob);
        else
            data.Y = data.Y;
        end
        % Check the W variable
        if size(data.W,3) > 1
            data.W = data.W(:,:,VoxelForThisJob);
        else
            data.W = data.W;
        end
        % Check the V variable
        if size(data.V,3) > 1
            data.V = data.V(:,:,VoxelForThisJob);
        else
            data.V = data.V;
        end
        % Check the covariates
        if size(data.COV,3) > 1
            data.COV = data.COV(:,VoxelForThisJob);
        else
            data.COV = data.COV;
        end
        
    case '74'
        if size(data.X,3) > 1
            data.X = data.X(:,:,VoxelForThisJob);
        else
            data.X = data.X;
        end
        % Check the M variable
        if size(data.M,3) > 1
            data.M = data.M(:,:,VoxelForThisJob);
        else
            data.M = data.M;
        end
        % Check the Y variable
        if size(data.Y,3) > 1
            data.Y = data.Y(:,:,VoxelForThisJob);
        else
            data.Y = data.Y;
        end
        % Check the covariates
        if size(data.COV,3) > 1
            data.COV = data.COV(:,VoxelForThisJob);
        else
            data.COV = data.COV;
        end
    case '75'
        if size(data.X,3) > 1
            data.X = data.X(:,:,VoxelForThisJob);
        else
            data.X = data.X;
        end
        % Check the M variable
        if size(data.M,3) > 1
            data.M = data.M(:,:,VoxelForThisJob);
        else
            data.M = data.M;
        end
        % Check the Y variable
        if size(data.Y,3) > 1
            data.Y = data.Y(:,:,VoxelForThisJob);
        else
            data.Y = data.Y;
        end
        % Check the covariates
        if size(data.COV,3) > 1
            data.COV = data.COV(:,VoxelForThisJob);
        else
            data.COV = data.COV;
        end
end

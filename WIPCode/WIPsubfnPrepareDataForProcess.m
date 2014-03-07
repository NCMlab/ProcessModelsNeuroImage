function WIPsubfnPrepareDataForProcess(ModelInfo)


BaseDir = ModelInfo.BaseDir;
Nvoxels = length(ModelInfo.Indices);
NJobSplit = ModelInfo.NJobSplit;
NSub = size(ModelInfo.data{1},1);

% how many jobs to split this into
NvoxelsPerJob = ceil(Nvoxels/NJobSplit);


% create a folder for data and results
c = clock;
OutFolderName = sprintf('%s_%s_%02d-%02d',ModelInfo.Tag,date,c(4),c(5))
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

% Create the bootstrap resamples
Samp = int16(zeros(NSub,Nboot));
for i = 1:Nboot;
    Samp(:,i) = WIPBootStrap(ModelInfo.STRAT,NSub);
end
ModelInfo.BootSamp = Samp;

AllParameters = cell(Nvoxels,1);

for i = 1:NJobSplit
    % This just makes sure that the loop is entered at least once
    if (NJobSplit > 1) & (i == NJobSplit)
        break
    end
    VoxelForThisJob = [(i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob];
    data = ModelInfo;
    for j = 1:size(data.data,2);
        switch size(data.data{j},2) 
            case Nvoxels
                data.data{j} = data.data{j}(:,VoxelForThisJob);
            case 1
            
            otherwise
                errordlg('Data the wrong size');
        end
    end
    data.Indices = VoxelForThisJob;
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
        fprintf(fid,'%s\n','addpath /share/studies/CogRes/Scripts/ProcessModelsNeuroImage/FinalCode');
%        fprintf(fid,'%s\n','addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode');
        fprintf(fid,'%s\n',['subfnVoxelWiseProcessBatch(''' InDataPath ''');']);
        fprintf(fid,'exit\n');
        fprintf(fid,'EOF\n');
        fclose(fid);
        %    Str = ['! qsub  ' jobPath];
        Str = ['! qsub -q short.q -p -10 -e ' JobOutputFolder ' -o ' JobOutputFolder ' -l mem_free=500M ' jobPath];
        unix(Str);
    else
        subfnVoxelWiseProcessBatch(InDataPath);
    end
end
if NJobSplit > 1
    % now run the last chunk of data
    VoxelForThisJob =[(i-1)*NvoxelsPerJob + 1:Nvoxels];
    data = AllData;
    data = ModelInfo;
    for j = 1:size(data.data,2);
        switch size(data.data{j},2) 
            case Nvoxels
                data.data{j} = data.data{j}(:,VoxelForThisJob);
            case 1
            
            otherwise
                errordlg('Data the wrong size');
        end
    end
    data.Indices = VoxelForThisJob;
    
    %Parameters = subfnVoxelWiseProcessBatch(temp,ModelNum,Nboot,Thresholds);
    InTag = sprintf('data_%04d',i);
    InDataPath = fullfile(DataFolder,InTag);
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
    fprintf(fid,'%s\n','addpath /share/studies/CogRes/Scripts/ProcessModelsNeuroImage/FinalCode');
    %fprintf(fid,'%s\n','addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode');
    fprintf(fid,'%s\n',['subfnVoxelWiseProcessBatch(''' InDataPath ''');']);
    fprintf(fid,'exit\n');
    fprintf(fid,'EOF\n');
    fclose(fid);
    Str = ['! qsub -q short.q -e ' JobOutputFolder ' -o ' JobOutputFolder ' -l mem_free=500M ' jobPath];
    unix(Str);
end

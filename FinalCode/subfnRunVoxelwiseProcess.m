function subfnRunVoxelwiseProcess(AllData,AnalysisParameters)
addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode

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
JobFolder = fullfile(OutFolder,'jobs');
if ~exist(OutFolder)
    mkdir(OutFolder);
end
% create a folder for all the cluster job files
if ~exist(JobFolder)
    mkdir(JobFolder);
end

% save the analysis parameters 
Str = ['save ' fullfile(OutFolder,'AnalysisParameters') ' AnalysisParameters'];
eval(Str);

AllParameters = cell(Nvoxels,1);

for i = 1:NJobSplit - 1
    
    VoxelForThisJob =[(i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob];
    data = AllData;
    data.Nboot = data.Nboot;
    switch ModelNum
        case '4'
            data.M = AllData.M(:,:,VoxelForThisJob);
        case '14'
            data.M = AllData.M(:,:,VoxelForThisJob);
            data.V = AllData.V(:,:,VoxelForThisJob);
    end
    data.Indices = AllData.Indices(VoxelForThisJob);
    %Parameters = subfnVoxelWiseProcessBatch(temp,ModelNum,Nboot,Thresholds);
    InTag = sprintf('data_%04d',i);
    InDataPath = fullfile(OutFolder,InTag);
    Str = ['save ' InDataPath ' data  '];
    eval(Str);
    jobPath = fullfile(OutFolder,sprintf('job_%04d.sh',i));
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
    Str = ['! qsub -q short.q -p -10 -e ' JobFolder ' -o ' JobFolder ' -l mem_free=500M -l h_vmem=1G ' jobPath];
    unix(Str);
end
% now run the last chunk of data
VoxelForThisJob =[(i)*NvoxelsPerJob + 1:Nvoxels];
data = AllData;
data.M = AllData.M(:,:,VoxelForThisJob);
% keep track of which voxels are being processed 
data.Indices = AllData.Indices(VoxelForThisJob);
%Parameters = subfnVoxelWiseProcessBatch(temp,ModelNum,Nboot,Thresholds);
InTag = sprintf('data_%04d',i+1);
InDataPath = fullfile(OutFolder,InTag);
Str = ['save ' InDataPath ' data '];
eval(Str);

jobPath = fullfile(OutFolder,sprintf('job_%04d.sh',i+1));
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
Str = ['! qsub -q short.q -e ' JobFolder ' -o ' JobFolder ' -l mem_free=500M -l h_vmem=1G ' jobPath];
unix(Str);
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


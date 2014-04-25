function WIPsubfnPrepareDataForProcessPERMUTE(ModelInfo)


BaseDir = ModelInfo.BaseDir;
Nvoxels = length(ModelInfo.Indices);
NJobSplit = ModelInfo.NJobSplit;
NSub = size(ModelInfo.data{1},1);
Nperm = ModelInfo.Nperm;
Nboot = ModelInfo.Nboot;


% create a folder for data and results
fprintf(1,'Creating the output folders.\n');
c = clock;
OutFolderName = sprintf('%s_%s_%02d-%02d',ModelInfo.Tag,date,c(4),c(5));
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
ModelInfo.Boot = Samp;

% NORMALIZE THE data before saving
fprintf(1,'Normalizing the data.\n');
NormtempData = ModelInfo.data;
for i = 1:ModelInfo.Nvar
    tempData = ModelInfo.data{i};
    [m n] = size(tempData);
    for j = 1:n
        tempData(:,j) = (tempData(:,j) - mean(tempData(:,j)))./std(tempData(:,j));
    end
    NormtempData{i} = tempData;
end
clear tempData
ModelInfo.data = NormtempData;
ModelInfo.OutFolder = OutFolder;

% Since this is the PERMUTATION the data only needs to be written once
% because it is not being broken into chunks
InTag = sprintf('AllData');
InDataPath = fullfile(DataFolder,InTag);
Str = ['save ' InDataPath ' ModelInfo  '];
eval(Str);
AnalysisParameters = ModelInfo;
AnalysisParameters.data = [];
% Save the analysis parameters which may be useful later
InParametersPath = fullfile(DataFolder,'AnalysisParameters');
Str = ['save ' InParametersPath ' AnalysisParameters '];
eval(Str);
% save the data and analysis parameters
InDataPath = fullfile(DataFolder,'ModelInfo');
Str = ['save ' InDataPath ' ModelInfo '];
eval(Str);

%% CALCULATE THE POINT ESTIMATE

  jobPath = fullfile(InJobFolder,sprintf('job_PointEstimate.sh',i));
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
        fprintf(fid,'%s\n',['WIPsubfnVoxelWiseProcessPERMUTE(''' InDataPath ''',''' num2str(0) ''',''' num2str(0) ''');']);
        fprintf(fid,'exit\n');
        fprintf(fid,'EOF\n');
        fclose(fid);
        %    Str = ['! qsub  ' jobPath];
        Str = ['! qsub -q short.q -p -10 -e ' JobOutputFolder ' -o ' JobOutputFolder ' -l mem_free=4G ' jobPath];
        unix(Str);
%% CALCULATE THE PERMUTATIONS
for i = 1:NJobSplit
    % If there is only ONE job specified then DO NOT send this to the
    % clusterwhos
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
        fprintf(fid,'%s\n',['WIPsubfnVoxelWiseProcessPERMUTE(''' InDataPath ''',''' num2str(i)  ''',''' num2str(Nperm/NJobSplit) ''');']);
        fprintf(fid,'exit\n');
        fprintf(fid,'EOF\n');
        fclose(fid);
        %    Str = ['! qsub  ' jobPath];
        Str = ['! qsub -q short.q -p -10 -e ' JobOutputFolder ' -o ' JobOutputFolder ' -l mem_free=4G ' jobPath];
        unix(Str);
    else
        subfnVoxelWiseProcessBatch(InDataPath);
    end
end



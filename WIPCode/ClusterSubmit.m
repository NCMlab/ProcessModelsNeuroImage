%% NOTES
% The program can split the data along voxels or along resamples. Splitting
% along voxels creates much smaller chunks of data and result files.
% Splitting along resamples creates data chunks of the entire brain which
% can be quite large and "clog" the memory of a computer. The problem is
% that it may be of interest to perform resampling the same way for every
% voxel. In that case the resamples can be determined before the data is
% split into chunks of voxels then this same resampling can be used across
% all data chunks.





%%
%BaseDir = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/Model75_XAgeGroup_MpGM_MfMRI_YPerf_COVNboot0_25-Feb-2014_11-46';
SelectedDir = spm_select(1,'dir');
Nboot = 5000;
NJobs = 1000;
NbootSplit = Nboot/NJobs;
% Calculate and save the point estimate results
cd(SelectedDir)
load AllData
% Calculating point estimate results before launching the boot strap
% process
OutFile = fullfile(SelectedDir,sprintf('PointEstimateResults.mat'));
if exist(OutFile,'file')
    fprintf(1,'Point estimates already calculated!\n');
else
    fprintf(1,'Calculating the point estimates\n');
    tic
    PointEstResults = WIPsubfnVoxelWiseProcessBatch(AllData);
    t = toc;
    fprintf(1,'Done in: %0.4f\n',t);
    fprintf(1,'Saving point estimates\n');
    SinglePointEstimate = PointEstResults{1}; 
    % Save the point estimates
    Str = sprintf('save %s PointEstResults',OutFile);
    eval(Str)
end

JobFolder = fullfile(SelectedDir,'Jobs');
JobOutputFolder = fullfile(SelectedDir,'JobOutput');
mkdir(JobFolder);
mkdir(JobOutputFolder);

for i = 1:NJobs
        jobPath = fullfile(JobFolder,sprintf('job_%04d.sh',i));
        fid = fopen(jobPath,'w');
        fprintf(fid,'#!/bin/bash\n');
        fprintf(fid,'#PBS -N Matlab\n');
        fprintf(fid,'#PBS -m be\n');
        fprintf(fid,'#PBS -j oe\n');
        fprintf(fid,'cd %s\n',SelectedDir);
        fprintf(fid,'/usr/local/matlab/bin/matlab -nodisplay << EOF\n');
        fprintf(fid,'%s\n','addpath /share/data/data5/spm8');
        fprintf(fid,'%s\n','addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode');
        fprintf(fid,'%s\n','addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/WIPCode');
        
        fprintf(fid,'%s\n',['WIPProcessBatch(''' SelectedDir ''',''' num2str(i) ''',''' num2str(NbootSplit) ''');']);
        fprintf(fid,'exit\n');
        fprintf(fid,'EOF\n');
        fclose(fid);
        %    Str = ['! qsub  ' jobPath];
        Str = ['! qsub -p -10 -e ' JobOutputFolder ' -o ' JobOutputFolder ' -l mem_free=500M ' jobPath];
        unix(Str);
end


function AddMorePermutations(ResultsFolder,TotalPerm,NPermPerJob)
% How many permutations have been done?
Results = fullfile(ResultsFolder,'Results');
Fres = dir(fullfile(ResultsFolder,'Results','Permute_count*.mat'));
NpermDone = length(Fres)
% What is the max permutation number
LastName = Fres(end).name;
FindUnder = strfind(LastName,'_')
MaxPerm = str2double(LastName(FindUnder(2)-4:FindUnder(2)-1));
% How many more are needed?
ToDoPerm = TotalPerm - NpermDone;

InJobFolder = fullfile(ResultsFolder,'JobFiles');
InDataPath = fullfile(ResultsFolder,'data','ModelInfo');
JobOutputFolder = fullfile(ResultsFolder,'JobOutput');

NJobSplit = ceil(ToDoPerm/NPermPerJob);

% Create a list of submitted jobs to be used for executing
% clean up commands after the jobs have finished.
WaitList = '';

% Submit the point estimate job
% Create the cluster submission job
jobPath = fullfile(InJobFolder,sprintf('PointEst_job.sh'));
fid = fopen(jobPath,'w');

% Create string of the command to be run with MatLab
Command = sprintf('VoxelWiseProcessPermute(''%s'',''%d'',''%d'');',InDataPath,0,0);
% Create the cluster submission script
CreateClusterJobFile(Command,fid)
% Submit the script to the cluster
JobName = SubmitClusterJob(jobPath,JobOutputFolder);
            

% Add job to wait list
WaitList = sprintf('%s,%s',WaitList,JobName);

% Submit all of the permutation jobs
for i = 1:NJobSplit
    % Create the cluster submission job
    jobPath = fullfile(InJobFolder,sprintf('job_%04d.sh',i));
    fid = fopen(jobPath,'w');
    % Create string of the command to be run with MatLab
    Command = sprintf('VoxelWiseProcessPermute(''%s'',''%d'',''%d'');',InDataPath,i,NPermPerJob);
    % Create the cluster submission script
    CreateClusterJobFile(Command,fid)
    % Submit the script to the cluster
    JobName = SubmitClusterJob(jobPath,JobOutputFolder);
    % Add job to wait list
    WaitList = sprintf('%s,%s',WaitList,JobName);
end
            
% Write out the resultant images
% Command = sprintf('WriteOutResults(''%s'')',OutFolder);
% jobPath = fullfile(InJobFolder,sprintf('WriteOutJob.sh'));
% fid = fopen(jobPath,'w');
% CreateClusterJobFile(Command,fid)
% % Submit the job with the wait list.
% SubmitClusterJob(jobPath,JobOutputFolder,WaitList);

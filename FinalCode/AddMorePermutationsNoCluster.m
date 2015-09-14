function AddMorePermutationsNoCluster(ResultsFolder,TotalPerm)
% How many permutations have been done?
Results = fullfile(ResultsFolder,'Results');
Fres = dir(fullfile(ResultsFolder,'Results','Path_count*.mat'));
NpermDone = length(Fres);
% What is the max permutation number
LastName = Fres(end).name;
Fdot = strfind(LastName,'_');
MaxPerm = str2double(LastName(end-7:end-4));
% How many more are needed?
ToDoPerm = TotalPerm - NpermDone;

%InJobFolder = fullfile(ResultsFolder,'JobFiles');
InDataPath = fullfile(ResultsFolder,'data','ModelInfo');
load(InDataPath)
%JobOutputFolder = fullfile(ResultsFolder,'JobOutput');


VoxelWiseProcessPermute(InDataPath,1,ToDoPerm,ModelInfo.TFCEparams,MaxPerm + 1);
            
% Write out the resultant images
% Command = sprintf('WriteOutResults(''%s'')',OutFolder);
% jobPath = fullfile(InJobFolder,sprintf('WriteOutJob.sh'));
% fid = fopen(jobPath,'w');
% CreateClusterJobFile(Command,fid)
% % Submit the job with the wait list.
% SubmitClusterJob(jobPath,JobOutputFolder,WaitList);

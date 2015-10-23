function AddMorePermutationsNoCluster(InFolder,TotalPerm)
% How many permutations have been done?
ResultsFolder = fullfile(InFolder,'Results');
Fres = dir(fullfile(ResultsFolder,'Permute_count*.mat'));
cd(ResultsFolder);
% From all of these files find out how many permutations have been
% completed
NpermDone  = 0;
for i = 1:length(Fres)
    temp = load(Fres(i).name);
    NpermDone = NpermDone + size(temp.MaxBeta,3);
end

ToDoPerm = TotalPerm - NpermDone;

%InJobFolder = fullfile(ResultsFolder,'JobFiles');
InDataPath = fullfile(InFolder,'data','ModelInfo');
load(InDataPath)
%JobOutputFolder = fullfile(ResultsFolder,'JobOutput');



ChunkSize = 2;
for i = 1:ceil(ToDoPerm/ChunkSize)
    VoxelWiseProcessPermute(InDataPath,ChunkSize,ModelInfo.TFCEparams,NpermDone+(i-1)*ChunkSize+1);
end


% Write out the resultant images
% Command = sprintf('WriteOutResults(''%s'')',OutFolder);
% jobPath = fullfile(InJobFolder,sprintf('WriteOutJob.sh'));
% fid = fopen(jobPath,'w');
% CreateClusterJobFile(Command,fid)
% % Submit the job with the wait list.
% SubmitClusterJob(jobPath,JobOutputFolder,WaitList);

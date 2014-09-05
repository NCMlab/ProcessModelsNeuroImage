function ResubmitJobs(InFolder)
% Sometimes jobs fail. This function checks to see which results are
% missing and resubmits the jobs. This of course assumes that the user has
% solved or fixed the problem that caused teh job to fail in the first
% place.
if nargin == 0
    InFolder = spm_select(1,'dir','Select folder');
end
% Load the ModelInfo
load(fullfile(InFolder,'data','ModelInfo'))
JobOutputFolder = fullfile(InFolder,'JobOutput');

ResultFiles = dir(fullfile(InFolder,'Results','Permute*.mat'));

NotFinishedJobs = ones(ModelInfo.NJobSplit,1);

NFiles = length(ResultFiles);
if ModelInfo.NJobSplit ~= NFiles
    % Not all jobs finished!
    for i = 1:NFiles
        
        Funder = findstr(ResultFiles(i).name,'_');
        Fdot = findstr(ResultFiles(i).name,'.');
        FileCount = str2double(ResultFiles(i).name(Funder(1)+6:Funder(2)-1));
        
        NotFinishedJobs(FileCount) = 0;
    end
end
FToDo = find(NotFinishedJobs);
for i = 1:length(FToDo)
    jobPath = fullfile(InFolder,'JobFiles',sprintf('job_%04d.sh',FToDo(i)));
    [JobName, ClusterCommand] = SubmitClusterJob(jobPath,JobOutputFolder);
end

fprintf(1,'%s\n',InFolder)
fprintf(1,'%d jobs resubmitted out of %d\n',length(FToDo),ModelInfo.NJobSplit)
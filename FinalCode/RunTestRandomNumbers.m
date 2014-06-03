JobOutputFolder='/home/js2746/DropBox/SteffenerColumbia/Projects/iLS/TestNum'

for i = 1:5
    jobPath = fullfile(JobOutputFolder,sprintf('job_%04d.sh',i));
    fid = fopen(jobPath,'w');
    Command = sprintf('TestRandomNumbersOnCluster');
    % Create the cluster submission script
    CreateClusterJobFile(Command,fid)
    % Submit the script to the cluster
    JobName = SubmitClusterJob(jobPath,JobOutputFolder);
end
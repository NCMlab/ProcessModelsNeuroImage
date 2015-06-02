function [JobName, ClusterCommand] = SubmitClusterJob(jobPath,JobOutputFolder,WaitList)
% This is the function that makes the qsub call for submitting a job to the
% cluster.
% Return the job name so it can be used with a wait command for executing
% cluster jobs that require other jobs to finish first.
% 
% This function needs some optimization. Specifically, I need to learn more
% about the memory management and use of the queues. It would be nice t
% estimate the memory requirements of the jobs based o n the data size and
% the number of jobs the analysis is split into. Such an approach would
% optimize the cluster use. 


if nargin == 2
    [PathName, JobName] = fileparts(jobPath);
    % Submit the job and explictly name the job.
    ClusterCommand = sprintf('! qsub -q all.q -p -50 -l h_vmem=8G,mem_free=6G,h_stack=256M -e %s -o %s -N %s %s',...
        JobOutputFolder,JobOutputFolder,JobName,jobPath);
   % Str = sprintf('! qsub -e %s -o %s -N %s %s',...
   %     JobOutputFolder,JobOutputFolder,JobName,jobPath);
    
    % Execute the job submission command
    unix(ClusterCommand);
elseif nargin == 3
    % A wait list was entered
    ClusterCommand = sprintf('! qsub -hold_jid %s -p -10 -e %s -o %s %s',...
        WaitList,JobOutputFolder,JobOutputFolder,jobPath);
   % Str = sprintf('! qsub -hold_jid %s -e %s -o %s %s',...
   %     WaitList,JobOutputFolder,JobOutputFolder,jobPath);

    unix(ClusterCommand);
end


%% THE FOLLOWING MAY LIMIT THE NUMBER OF CONCURRENT JOBS
% qsub script.sh -t 0-99%10
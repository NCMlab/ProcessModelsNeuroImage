function SubmitClusterJob(jobPath,JobOutputFolder)
% This is the function that makes the qsub call for submitting a job to the
% cluster. 

Str = ['! qsub -q short.q -p -10 -e ' JobOutputFolder ' -o ' JobOutputFolder ' -l mem_free=500M ' jobPath]
unix(Str);
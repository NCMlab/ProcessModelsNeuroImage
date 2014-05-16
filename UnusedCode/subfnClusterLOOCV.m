function subfnClusterLOOCV(dataFile, NSub, OutDir, NPCs)


% Create the job file
jobPath = fullfile(OutDir,sprintf('job_LOOCV.sh'));
fid = fopen(jobPath,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'#PBS -N Matlab\n');
fprintf(fid,'#PBS -m be\n');
fprintf(fid,'#PBS -j oe\n');
fprintf(fid,'cd %s\n',OutDir);
fprintf(fid,'/usr/local/matlab/bin/matlab -nodisplay << EOF\n');
fprintf(fid,'%s\n','addpath /share/data/data5/spm8');
fprintf(fid,'%s\n','addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode');
fprintf(fid,'%s\n',['subfnCalculateLOOCV(''' dataFile ''',''' num2str(NPCs) ''');']);
fprintf(fid,'exit\n');
fprintf(fid,'EOF\n');
fclose(fid);

% This -t command for qsub submits the jobs as an array. The only
% difference is that only one job file is created. This uses an
% environmental variable to determin which call is being made andthat
% becomes the left out subject.
Str = ['! qsub -q short.q -p -10 -e ' OutDir ' -o ' OutDir ' -l mem_free=500M -t 1-' num2str(NSub) ' ' jobPath];
fprintf(1,'%s\n',Str);
unix(Str);
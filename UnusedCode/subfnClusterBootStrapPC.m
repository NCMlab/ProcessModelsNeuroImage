function subfnClusterBootStrapPC(data,Nboot,selected_PCs,NTimes,OutDir)
% Create the data file
dataFile = fullfile(OutDir,[sprintf('BootData')]);
Str = sprintf('save %s %s %s %s %s %s',dataFile,'data', 'Nboot','selected_PCs','OutDir');
eval(Str);

    % Create the job file
    jobPath = fullfile(OutDir,sprintf('BOOTStrap_job.sh'));
    fid = fopen(jobPath,'w');
    fprintf(fid,'#!/bin/bash\n');
    fprintf(fid,'#PBS -N Matlab\n');
    fprintf(fid,'#PBS -m be\n');
    fprintf(fid,'#PBS -j oe\n');
    fprintf(fid,'cd %s\n',OutDir);
    fprintf(fid,'/usr/local/matlab/bin/matlab -nodisplay << EOF\n');
    fprintf(fid,'%s\n','addpath /share/data/data5/spm8');
    fprintf(fid,'%s\n','addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode');
    fprintf(fid,'%s\n',['subfnCalculateBootStrapPC(''' dataFile ''');']);
    fprintf(fid,'exit\n');
    fprintf(fid,'EOF\n');
    fclose(fid);
    %    Str = ['! qsub  ' jobPath];
    Str = ['! qsub -q short.q -p -10 -e ' OutDir ' -o ' OutDir ' -l mem_free=500M -t 1-' num2str(NTimes) ' ' jobPath];
    unix(Str);

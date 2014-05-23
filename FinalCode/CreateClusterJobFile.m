function CreateClusterJobFile(Command,fid)
% This is a general use function for creating a shell script that will be
% submitted to a cluster. The only input is the actual MatLab command and
% the inputs to that command. The inputs are strings and NOT variables
% containing data.
% The purpose of this file is so that if the location of the MatLab
% install, or SPM, changes then only this file needs to be updated.
% 
% Here is an example command:
% 'subfnVoxelWiseProcessBatch(''' InDataPath ''');'
% It is a string and the information passed to it is a string in quotes. 
% For this case it is the path to a data file that will be used by the
% program.
        
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'#PBS -N Matlab\n');
fprintf(fid,'#PBS -m be\n');
fprintf(fid,'#PBS -j oe\n');
%fprintf(fid,'cd %s\n',BaseDir);
% Start matlab using the nodisplay option
fprintf(fid,'/usr/local/MATLAB/R2013b/bin/matlab -nodisplay << EOF\n');
% The pathe to the spm software
fprintf(fid,'%s\n','addpath /usr/local/SPM/v8');
% The path to the Process toolbox
fprintf(fid,'%s\n','addpath /home/js2746/DropBox/SteffenerColumbia/Scripts/ProcessModelsNeuroImage/FinalCode');
fprintf(fid,'%s\n',Command);
fprintf(fid,'exit\n');
fprintf(fid,'EOF\n');
fclose(fid);

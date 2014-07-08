% Write out basic regression maps using the minimum statistics correction
% for multiple comparisons from the permutation results.
%
% Add the ability to perform FDR correction for multiple comparisons when
% using bootstrapping. This requires the BCa Z and p maps to be saved.
%
% Add the ability to test multiple paths. 
%
% Include test for total effects with the mediation. This will need the
% addition of a third dimension to the DIRECT array.
% 
% Add the ability to test multiple paths in the model.
%
% Add a check for whether the data is "good" or not. This would be an
% implicit maskign sort of approach. At the simplest it should check for
% voxels with NaN and check whether there are "bad" subjects or bad
% "voxels."  Then the user may be provided with an option of removing the
% voxel or subject. If a subject/voxel with NaN values is to be retained
% then there needs to be a procedure for interpolating it's value.
%
% Have permutations saved to a file in small chunks instead of the cluster
% split size. This is valuable for processes that get interrupted. 
% This is now a necessity. The memory demands for multiple permutations is
% HUGE. So it is best to write each permutation to a file after completion.

% There are some careful interplays between the use of qsub and MatLab
% multithreading nature. By default MatLab will multithread. This puts it
% in competition (I think) with qsub which determines how to allocate
% resources. Qsub allocates resources based on the quantity of available
% memory.
% I have the question of whether the multithreading of MatLab and qsub may
% be competing and that qsub may overallocte resources to a single compute
% host.

% Using 8 threads a process took 224 seconds on the head node. I want to
% see if using 1 thread takes 1/8 of the time or not. 
% Using one thread took 240 seconds!

% I am having lots of problems committing jobs to the cluster due to memory
% management. Jobs get almost done and then are killed. I am guessing that
% this is due to the job exceeding the memory limits. It would be useful to
% estimate the required memory of a job and use that value in the
% submission of a qsub. The problem seems to be that the results structure
% increases in size with every voxel. Therefore, by the end of a process
% this structure can get TOO big. So if a single voxel is analyzed for the
% selected model then the expected memory requirements would be estimated
% from the number of voxels assigned to a job. Then recommendations to the
% user can be provided based on these estimates.


% Add robust regression option

function AddMorePermutationsNoClusterResourceCheck(InFolder,TotalPerm)
% Only run if there are adequate resources for MatLab

CPUlimit = 50;
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
ErrorCount = 0;
for i = 1:ceil(ToDoPerm/ChunkSize)
    while ErrorCount < 100
        try
            % Some times a permutation will not work and I don't know why. This
            % just catches the error and tries again.
            
            ThisChunkRunFlag = 0;
            while ~ThisChunkRunFlag
                [Status, RES] = unix('top -bn 4 -d 0.01 | grep ''^%Cpu'' | tail -n 1 | gawk ''{print $2+$4+$6}''');
                % Check to see if the results can be converted to anumber? If so then
                % it worked.
                if ~isempty(str2num(RES))
                    % Check to see if the CPU usage limit is reached
                    if str2num(RES) < CPUlimit
                        fprintf(1,'Within CPU (%0.1f) limits, running a permutation\n',str2num(RES));
                        ThisChunkRunFlag = 1;
                    else
                        fprintf(1,'Not enough resources (%0.1f > %d) at %s, checking again in 10 min.\n',str2num(RES),CPUlimit,datestr(now));
                        % wait for ten minutes and check again
                        pause(10*60)
                    end
                else
                    error('The CPU status check did not work. Please use AddMorePermutationsNoCluster.m');
                end
            end
            VoxelWiseProcessPermute(InDataPath,ChunkSize,ModelInfo.TFCEparams,NpermDone+(i-1)*ChunkSize+1);
            % This will end the while loop
            ErrorCount = 1000;
        catch
            ErrorCount = ErrorCount + 1;
            fprintf('There was an error, trying again.\n');
        end
    end
end

% Write out the resultant images
% Command = sprintf('WriteOutResults(''%s'')',OutFolder);
% jobPath = fullfile(InJobFolder,sprintf('WriteOutJob.sh'));
% fid = fopen(jobPath,'w');
% CreateClusterJobFile(Command,fid)
% % Submit the job with the wait list.
% SubmitClusterJob(jobPath,JobOutputFolder,WaitList);

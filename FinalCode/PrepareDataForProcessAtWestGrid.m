function OutFolder = PrepareDataForProcessAtWestGrid(ModelInfo)


% The folder for housing these analyses uses the provided Tag name along
% with the current date and time.
c = clock;

% Create the output folder name
OutFolderName = sprintf('%s_%s_%02d-%02d',ModelInfo.Tag,date,c(4),c(5));

% Create the full output folder name
OutFolder = fullfile(ModelInfo.BaseDir,OutFolderName);

% Actually create the folders if they do not exist already.
if ~exist(OutFolder)
    mkdir(OutFolder);
end

% Create a folder to store the data in
DataFolder = fullfile(OutFolder,'data');
if ~exist(DataFolder)
    mkdir(DataFolder);
end

% Is this a bootstrap analysis or a permutation analysis?
if ModelInfo.Nboot > 0 && ModelInfo.Nperm > 0
    errordlg('The model can apply bootstrapping OR permutation testing NOT both.');
elseif ModelInfo.Nboot > 0
    ModelType = 'bootstrap';
elseif ModelInfo.Nperm > 0
    ModelType = 'permutation';
end

% Is this analysis using the cluster?
%% Is this analysis using the cluster?
if ModelInfo.NJobSplit > 1
    % Since this analysis will use the cluster then setup the folders to
    % hold the job files and their output
    % Create the folder to hold the output from the jobs executed on the
    % cluster. These are especially useful for investigation of errors.
    JobOutputFolder = fullfile(OutFolder,'JobOutput');
    if ~exist(JobOutputFolder,'dir')
        mkdir(JobOutputFolder);
    end
    
    % Create the folder to hold the job files that get sent to the cluster
    % itself. These are useful if any jobs
    InJobFolder = fullfile(OutFolder,'JobFiles');
    if ~exist(InJobFolder,'dir')
        mkdir(InJobFolder);
    end
    
    switch ModelType
        case 'permute'
        case 'bootstrap'
            % For some reason the cd is causing problems!
            cd(JobOutputFolder)
            % Here the data is split and saved as a series of small files.
            
            
            % How many voxels will be in each split?
            % Note that the ceil function rounds UP. So after hundreds of
            % jobs there may not NEED to be 500 jobs but only 498.
            NvoxelsPerJob = ceil(ModelInfo.Nvoxels/ModelInfo.NJobSplit);
            
            % Cycle over the number of splits and submit each data chunk to
            % the cluster to be processed.
            % Create a list of submitted jobs to be used for executing
            % clean up commands after the jobs have finished.
            WaitList = '';
            % This counter is used because sometimes the final job(s) does
            % not have any voxels in it vbecause of the ceil round
            % performed just above.
            ActualNJobSplit = 0;
            for i = 1:ModelInfo.NJobSplit
                
                VoxelsForThisJob = [(i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob];
                
                
                % remove any voxels in this split which are beyond the
                % total number of voxels in the dataset. This is only a
                % concern for the last data split which may be a partial
                % smaller data set.
                VoxelsForThisJob(find(VoxelsForThisJob > ModelInfo.Nvoxels)) = [];
                
                % Check to see if this list is empty
                if ~isempty(VoxelsForThisJob)
                    ActualNJobSplit = ActualNJobSplit + 1;
                    % Create a version of the full data/model structure only containing the
                    % subset of voxels to be analyzed.
                    subModelInfo = ModelInfo;
                    subModelInfo.Indices = VoxelsForThisJob;
                    for j = 1:ModelInfo.Nvar
                        % Check to see which variables are multi-voxel variables and
                        % extract the subset of data
                        if size(subModelInfo.data{j},2) > 1
                            subModelInfo.data{j} = ModelInfo.data{j}(:,VoxelsForThisJob);
                        end
                    end
                    % Save this subset of data to a file
                    InTag = sprintf('data_%04d',i);
                    InDataPath = fullfile(DataFolder,InTag);
                    Str = ['save ' InDataPath ' subModelInfo  '];
                    eval(Str);
                end
            end
            % Write out model info
            ModelInfo.NJobSplit = ActualNJobSplit;
            Str = sprintf('save %s ModelInfo',fullfile(DataFolder,'ModelInfo'));
            eval(Str);
            
    end
end


% 
% %% 
% % Then the data needs to be transferred



% scp -p -P 1234 -r Model1_ACC_sRet_HighInter_5000boot_10-Jul-2014_08-49/ steffejr@127.0.0.1:~/Projects/iLS/.
% 
% % At this point the python script needs to be called on westgrid
% 
% 
% 
% 
% [JobName,ClusterCommand] = SubmitClusterJob(jobPath,JobOutputFolder);
% fprintf(1,'%s\n',ClusterCommand);
% 
% % Replace this with the tools from FieldTrip
% %[JobName, puttime] = qsubfeval('CycleOverVoxelsProcessBootstrap',InDataPath,'memreq',1024^3,...
% %    'timreq',24*60*60,'matlabcmd','/usr/local/MATLAB/R2013b/bin/matlab',...
% %    'jvm','no','queue','veryshort.q');
% 
% 
% % Add job to wait list
% WaitList = sprintf('%s,%s',WaitList,JobName);
% end
% 
% % Remove initial comma from the WaitList
% WaitList = WaitList(2:end);
% % The following command is submitted to the cluster and given a
% % list of all the previous processing jobs. This forces the
% % write out of the resultant images to
% % wait for all of the other jobs to finish first. These
% % clean up commands and writing out the images do not need to be
% % run on the cluster.
% % The advantage of using the cluster queueing is that this job
% % waits until other jobs finish. By having the cluster wait
% % it frees up the MatLab process on the head node.
% %
% % Write out the resultant images
% Command = sprintf('WriteOutResults(''%s'')',OutFolder);
% jobPath = fullfile(InJobFolder,sprintf('WriteOutJob.sh'));
% fid = fopen(jobPath,'w');
% CreateClusterJobFile(Command,fid)
% % Submit the job with the wait list.
% SubmitClusterJob(jobPath,JobOutputFolder,WaitList);
% 
% 
% % Save the data
% % This is actually redundent because the data is being split
% % and then the subsets saved. But it is very useful for
% % rerunning a job(s).
% ModelInfo.NJobSplit = ActualNJobSplit;
% InDataPath = fullfile(DataFolder,'ModelInfo');
% Str = ['save ' InDataPath ' ModelInfo '];
% eval(Str);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

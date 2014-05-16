function OutFolder = PrepareDataForProcess(ModelInfo)
% In order for this program to work in both a cluster environment and a
% single computer environment the data is saved and the analysis programs
% are given a path name. The data specified by the path is then loaded, the
% data analyzed and the results saved. WHen using a cluster environment the
% data can be split into chunks and multiple function calls are sent to teh
% cluster with each provided a path to a different chunk of data. Instead
% of passing data to each function call on the cluster, only the paths are
% passed.
% 
% TO DO
% The programs are also smart enough to know that if actual data is passed,
% instead of a parth, then the data is considered pre-loaded and processed.


% This function checks the specified model to ensure it is correctly
% specified. It also sets up the output folder structure for housing
% results. This is especially important for working with a cluster
% environment.
% 

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
else
    % This analysis is run on the host computer
    switch ModelType
        case 'bootstrap'
            % Need to split the data and cycle over voxels
            Results = CycleOverVoxelsProcessBootstrap(ModelInfo);
        case 'permutation'
            % The permutation analysis requires calculation of the point
            % estimate and then calculation of all of the permutations. The
            % point estimate is a full set of results and the permutations
            % only require a smaller subset to be calculated. Therefore,
            % the permutation test requires the saving of data to file.
            % The permutation test is also used to calculate the multiple
            % comparison corrected results across voxels. Therefore, it
            % does not make sense to run this type of analysis on a single
            % point data set.
            % 
            % Save the data 
            InDataPath = fullfile(DataFolder,'ModelInfo');
            Str = ['save ' InDataPath ' ModelInfo '];
            eval(Str);
            
            % Perform the point estimate calculation
            VoxelWiseProcessPermute(InDataPath,0,0)
            % Perform the permutation tests
            VoxelWiseProcessPermute(InDataPath,1,ModelInfo.Nperm)
            
    end
end
    
    


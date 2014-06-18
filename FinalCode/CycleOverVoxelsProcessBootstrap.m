function Results = CycleOverVoxelsProcessBootstrap(ModelInfo)
% This function cycles over all voxels in the data chunk and processes
% them. This data "chunk" may in fact be the entire data set.
% In order to maximize the re-use of code this function checks to see if it
% is passed a structure or a filename. If it is passed a structure the data
% in that structure is processed. If it is passed a file name, the data is
% loaded from the file and then processed. 
%

ClusterJobFlag = 0;
if ~isstruct(ModelInfo)
    % If a structure is NOT passed
    if exist([ModelInfo '.mat'],'file') == 2
        ClusterJobFlag = 1;
        % If this function is being used by a cluster and reading the input
        % data from a file, the output also needs to be saved to a file.
        % Create the Results folder for the data if it has not already been
        % created.
        [PathName, FileName] = fileparts(ModelInfo);
        [CurrentAnalysisDir] = fileparts(PathName);
        ResultsDir = fullfile(CurrentAnalysisDir,'Results');
        % Find out which data chunk this is
        FindUnder = findstr(FileName,'_');
        DataChunk = FileName(FindUnder+1:end);
        % If the results folder does NOT exist then create it
        if ~exist(ResultsDir,'dir')
            mkdir(ResultsDir)
        end

        M = load(ModelInfo);
        % Once the data is loaded rename it
        F = fieldnames(M);
        ModelInfo = getfield(M,F{1});
    else
        errordlg('Wrong data type passed.')
    end
end
   

% The number of voxels variable refers to the number of voxels being
% analyzed in this chunk. If the data is being split into chunks for
% analysis on a cluster then this number will differe from teh total number
% of voxels in the analysis.
Nvoxels = length(ModelInfo.Indices);
% Prepare the output structure
Results = cell(Nvoxels,1);
if Nvoxels > 1
    for i = 1:Nvoxels
        t = tic;
        
        % Extract the data for this voxel
        OneVoxelModel = ExtractDataFromVoxel(ModelInfo,i);
        % Perform the analysis for this voxel.
        % SOmetines the bootstrapping produces "bad" results where the BCa
        % CI cannot be calculated. Therefore, a bad voxel is reprocessed
        % again to see if it wotks better the second time. The while loop
        % will repeat retying up to 20 times. If it cannot be fixed it is
        % skipped.
        
        ResultsFlag = 1;
        count = 0;
        while ResultsFlag
            Results{i} = OneVoxelProcessBootstrap(OneVoxelModel);
            count = count + 1;
            if isempty(Results{i})
                fprintf(1,'Did not work!!!\n')
            elseif count == 20
                ResultsFlag = 0;
            else
                ResultsFlag = 0;
            end
        end
        fprintf(1,'%d of %d voxels in %0.2f seconds.\n',i,Nvoxels,toc(t));
    end
else
    % Extract the data for this voxel
    OneVoxelModel = ExtractDataFromVoxel(ModelInfo,1);
    % Perform the analysis for this voxel
    ResultsFlag = 1;
    count = 0;
    while ResultsFlag
        Results{1} = OneVoxelProcessBootstrap(OneVoxelModel);
        count = count + 1;
        if isempty(Results{1})
            fprintf(1,'Did not work!!!\n')
        elseif count == 20
            ResultsFlag = 0;
        else
            ResultsFlag = 0;
        end
    end 
end

% Once everything is done, return the results to the calling function
% if a structure of data was used as input. If the input was a path to a
% file of data then save the results to a file.
if ClusterJobFlag
    ResultsFile = fullfile(ResultsDir,sprintf('BootStrapResults_%s',DataChunk));
    Str = sprintf('save %s Results',ResultsFile);
    eval(Str)
end


whos
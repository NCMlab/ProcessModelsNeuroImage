function VoxelWiseProcessPermute(InDataPath,count,Nperm)
% Make sure the random number seed is reset for every function call. This
% avoids each cluster node starting with the same seed and producing the
% same results. An alternative futur direction is a precalculation of the permutations for
% storage in the ModelInfo structure. 
RandStream.setDefaultStream(RandStream('mt19937ar','Seed',sum(100*clock)));

fprintf(1,'Started at: %s\n',datestr(now));
tic

% load data
if nargin == 3
    % Load the data and other variables
    load(InDataPath)
    count = str2num(count);
    Nperm = str2num(Nperm);
elseif nargin == 2
    % If the number of permutations field is left blank then assume that
    % only point estimate is being calculated.
    load(InDataPath)
    count = str2num(count);
    Nperm = 0;
elseif nargin == 1
    % point estimate
    load(InDataPath)
    count = 0;
    Nperm = 0;
else
    error('Expected at least one input.');
end
% Ensure that the InDataPath actually contained the ModelInfo structure
if ~exist('ModelInfo','var')
    errordlg('The data passed does not contain the ModelInfo structure.');
end

% extract the number of subjects
NSub = ModelInfo.NSub;

% Determine the sample order. This is where the permutation occurs.
if Nperm == 0
    % This is the point estimate
    Samp = [1:NSub]';
else
    % all permutations for this function call are specified here as
    % seperate columns.
    Samp = zeros(NSub,Nperm);
    for i = 1:Nperm
        Samp(:,i) = randperm(NSub)';
    end
end

% How many voxels 
Nvoxels = length(ModelInfo.Indices);
%%
if Nperm > 0
    
    % Create the structure to hold the permutation resample results.
    % Test one voxel to determine the correct size for all of the results.
    % The importance of this is because the user can specify any number of
    % paths to test. SOme of these paths may be simple, have no moderating
    % variables, while some may be more complex having interaction terms
    % which need to be probed.
    tempData = ModelInfo;
    
    % Since the data may contain multiple voxels for any of the nodes in
    % the path model this extracts only a single voxel and runs the process
    % modeling on it.
    tempData.data = zeros(ModelInfo.NSub,ModelInfo.Nvar);
    
    % Cycle over each data variable and pulls out the first column whether
    % it is multi-column of only has a single column.
    % extract the data
    for j = 1:ModelInfo.Nvar
        % CHeck to see if this variable has more then one column
        if size(ModelInfo.data{j},2) > 1
            tempData.data(:,j) = ModelInfo.data{j}(:,1);
        else
            tempData.data(:,j) = ModelInfo.data{j};
        end
    end
    % Fit the model for this test data
    TestResults = FitProcessModel(tempData);
    
    % Once the test results are calculated a structure is created to hold
    % the path results.
    PermResults = {};
    % Here any number of fields can be extracted by adding their names to
    % this list of Field Names. This would allow permutation tests on any
    % of the other results. Our main interest however is to apply
    % permutation testing to just the Paths.
    FieldNames = {'Paths'};
    for k = 1:length(FieldNames)
        Value = getfield(TestResults,FieldNames{k});
        if iscell(Value)
            BlankValue = cell(size(Value,1),Nvoxels);
        else
            BlankValue = zeros([size(Value) Nvoxels]);
        end
        PermResults = setfield(PermResults,FieldNames{k},BlankValue);
    end
    % determine the size of the matrix of regression parameters. The reason
    % this is not simply the number of variables entered is because there
    % can be any number of interactions in the model.
    [Mbeta, Nbeta] = size(TestResults.B);
    
    % Create the empty array to hold the parameter estimates
    PermResults.beta = zeros(Mbeta, Nbeta, Nvoxels, Nperm);
    
    % Determine the size of the calculated path values
    [PathSize1 PathSize2] = size(TestResults.Paths{1});
    
    % Determine the number of paths
    NPaths = length(TestResults.Paths);
    
    % THe permutation analysis uses the minimum statistic method for
    % correcting for multiple comparisons. This approach creates a matrix
    % of maximum and minimum values for each permutation.
    %
    % Nichols TE, Holmes AP. Nonparametric permutation tests for functional neuroimaging: 
    % a primer with examples.  Hum Brain Mapp. Department of Biostatistics, 
    % University of Michigan, Ann Arbor, Michigan; Robertson Centre for Biostatistics, 
    % Department of Statistics, University of Glasgow, Scotland, United Kingdom; 
    % Wellcome Department of Cognitive Neurology, Institute of Neurology, London, 
    % United Kingdom; 2002 Jan;15(1):1?25. 
    
    MaxB = zeros([size(TestResults.B) Nperm]);
    MinB = zeros([size(TestResults.B) Nperm]);
    MaxPaths = zeros(PathSize1,PathSize2,NPaths,Nperm);
    MinPaths = zeros(PathSize1,PathSize2,NPaths,Nperm);
end
%%

% How long did it take to prepare the data
fprintf(1,'Data prepared in %0.2f s.\n',toc);

% Feedback for the user
fprintf(1,'Starting first permutation at: %s\n',datestr(now));

% Cycle over 
for k = 1:size(Samp,2)
    tic
    for i = 1:Nvoxels
        % for this voxel
        tempData = ModelInfo;
        tempData.data = zeros(ModelInfo.NSub,ModelInfo.Nvar);
        % extract the data
        for j = 1:ModelInfo.Nvar
            if size(ModelInfo.data{j},2) > 1
                tempData.data(:,j) = ModelInfo.data{j}(:,i);
            else
                tempData.data(:,j) = ModelInfo.data{j};
            end
        end
        % Shuffle the first column
        tempData.data(:,1) = tempData.data(Samp(:,k),1);
        Results = WIPsubfnFitModel(tempData);
        if Nperm == 0
            % this is the point estimate, keep everything
            Parameters{i} = Results;
        else
            % only keep the path estimates
            PermResults.Paths{i} = Results.Paths;
            PermResults.beta(:,:,i) = Results.beta;
        end
        
        % Doing this is real nice to see progress but it SLOWS it down A LOT!!!
        % if ~mod(i,1000)
        %     fprintf(1,'Finished voxel %d of %d in %0.2f s.\n',i,Nvox,toc);
        % end
    end
    if Nperm > 0
        fprintf(1,'\tFinding the maxima and minima from this permutaion.\n');
        clear Parameters
        % Find the max and min values from the permutation
        for i = 1:NPaths
            ThisPath = zeros(PathSize1,PathSize2,Nvoxels);
            for j = 1:Nvoxels
                ThisPath(:,:,j) = PermResults.Paths{j}{i};
            end
            for j = 1:PathSize1
                for m = 1:PathSize2
                    MaxPaths(j,m,i,k) = max(squeeze(ThisPath(j,m,:)));
                    MinPaths(j,m,i,k) = min(squeeze(ThisPath(j,m,:)));
                end
            end
        end
        % find the max and min beta values also
        for i = 1:Mbeta
            for j = 1:Nbeta
                MaxB(i,j,k) = max(squeeze(PermResults.beta(i,j,:)));
                MinB(i,j,k) = min(squeeze(PermResults.beta(i,j,:)));
            end
        end
    end
    fprintf(1,'Finished permutation %d of %d in %0.2f s.\n',k,Nperm,toc);
end
fprintf(1,'Saving data to file now.\n');
% Save results to file
if Nperm > 0
    [PathName FileName] = fileparts(InDataPath);
    [PathName FileName] = fileparts(PathName);
    ResultsFolder = fullfile(PathName,'Results');
    if ~exist(ResultsFolder,'dir')
        mkdir(ResultsFolder)
    end
    OutFile = fullfile(ResultsFolder,sprintf('Permute_count%04d_%dSamp',count,Nperm));
    Str = sprintf('save %s MaxPaths MinPaths MaxB MinB',OutFile);
    eval(Str)
else
    [PathName FileName] = fileparts(InDataPath);
    [PathName FileName] = fileparts(PathName);
    ResultsFolder = fullfile(PathName,'Results');
    if ~exist(ResultsFolder,'dir')
        mkdir(ResultsFolder)
    end
    OutFile = fullfile(ResultsFolder,'PointEstimate');
    Str = sprintf('save %s Parameters',OutFile);
    eval(Str)
end
    
    
    
    
    
    
    
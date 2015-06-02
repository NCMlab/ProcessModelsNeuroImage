function VoxelWiseProcessPermute(InDataPath,count,Nperm,TFCEparams)
% Make sure the random number seed is reset for every function call. This
% avoids each cluster node starting with the same seed and producing the
% same results. An alternative futur direction is a precalculation of the permutations for
% storage in the ModelInfo structure. 
%rng('shuffle','multFibonacci');

setenv('FSLOUTPUTTYPE','NIFTI');
path1 = getenv('PATH');
if isempty(strfind(path1,'fsl'))
    path1 = [path1 ':/usr/local/fsl/bin'];
    setenv('PATH', path1);
end

%RandStream('mt19937ar','Seed',sum(100*clock));
rng shuffle

fprintf(1,'Started at: %s\n',datestr(now));
tic
fprintf(1,'%s\n',InDataPath);
% load data
if nargin == 4
    % Load the data and other variables
    load(InDataPath)
    % The cluster function calls pass strings
    if isstr(count); count = str2num(count); end
    if isstr(Nperm); Nperm = str2num(Nperm); end
elseif nargin == 3
    % Load the data and other variables
    load(InDataPath)
    % The cluster function calls pass strings
    if isstr(count); count = str2num(count); end
    if isstr(Nperm); Nperm = str2num(Nperm); end
    TFCEparams = [2 0.5 6];
elseif nargin == 2
    % If the number of permutations field is left blank then assume that
    % only point estimate is being calculated.
    load(InDataPath)
    if isstr(count); count = str2num(count); end
    Nperm = 0;
     TFCEparams = [2 0.5 6];
elseif nargin == 1
    % point estimate
    load(InDataPath)
    count = 0;
    Nperm = 0;
     TFCEparams = [2 0.5 6];
else 
    error('Expected at least one input.');
end


% Ensure that the InDataPath actually contained the ModelInfo structure
if ~exist('ModelInfo','var')
    errordlg('The data passed does not contain the ModelInfo structure.');
end

% extract the number of subjects
Nsub = ModelInfo.Nsub;

% Determine the sample order. This is where the permutation occurs.
if Nperm == 0
    % This is the point estimate
    Samp = [1:Nsub]';
else
    % all permutations for this function call are specified here as
    % seperate columns.
    Samp = zeros(Nsub,Nperm);
    for i = 1:Nperm
        Samp(:,i) = randperm(Nsub)';
    end
end

% How many voxels 
Nvoxels = length(ModelInfo.Indices);
%%
[PathName FileName] = fileparts(InDataPath);
[PathName FileName] = fileparts(PathName);
ResultsFolder = fullfile(PathName,'Results');
if ~exist(ResultsFolder,'dir')
    mkdir(ResultsFolder)
end


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
    tempData.data = zeros(ModelInfo.Nsub,ModelInfo.Nvar);
    
    % Cycle over each data variable and pulls out the first column whether
    % it is multi-column of only has a single column.
    % extract the data
    for j = 1:ModelInfo.Nvar
        % Check to see if this variable has more then one column
        if size(ModelInfo.data{j},2) > 1
            tempData.data(:,j) = ModelInfo.data{j}(:,j);
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
    
    clear BlankValue;
    % determine the size of the matrix of regression parameters. The reason
    % this is not simply the number of variables entered is because there
    % can be any number of interactions in the model.
    [Mbeta, Nbeta] = size(TestResults.B);
    
    % Create the empty array to hold the parameter estimates
    PermResults.beta = zeros(Mbeta, Nbeta, Nvoxels);
    PermResults.t = zeros(Mbeta, Nbeta, Nvoxels);
    % Determine the size of the calculated path values
    [PathSize1, PathSize2] = size(TestResults.Paths{1});
    
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
    
    MaxBeta = zeros([size(TestResults.B) Nperm]);
    MinBeta = zeros([size(TestResults.B) Nperm]);
    MaxPaths = zeros(PathSize1,PathSize2,NPaths,Nperm);
    MinPaths = zeros(PathSize1,PathSize2,NPaths,Nperm);
    TFCEtMax = zeros([size(TestResults.B) Nperm]);
    TFCEtMin = zeros([size(TestResults.B) Nperm]);
    TFCEpathsMax = zeros(PathSize1,PathSize2,NPaths,Nperm);
    TFCEpathsMin = zeros(PathSize1,PathSize2,NPaths,Nperm);
end
%%
tempData = ModelInfo;  
% How long did it take to prepare the data
fprintf(1,'Data prepared in %0.2f s.\n',toc);

% Feedback for the user
fprintf(1,'Starting first permutation at: %s\n',datestr(now));

% Create the structures for saving the temp files for use with the TFCE
% program in FSL
Vo = tempData.DataHeader;
Vo.fname = fullfile(tempData.BaseDir,'TestTemp.nii');
Vo.dt = [64 0];
tempI = zeros(Vo.dim);
% Create a temp mask image
Im = zeros(Vo.dim);
Im(tempData.Indices) = 1;
Vm = Vo;
Vm.fname = fullfile(tempData.BaseDir,'tempMask.nii');
spm_write_vol(Vm,Im);
outFile = fullfile(tempData.BaseDir,'tempTFCEout.nii');

% Cycle over 
for k = 1:size(Samp,2)
    tic
    % Cycle over all voxels
    for i = 1:Nvoxels
        % for this voxel
        % reset the data part
        tempData.data = zeros(ModelInfo.Nsub,ModelInfo.Nvar);
        % extract the data
        for j = 1:ModelInfo.Nvar
            if size(ModelInfo.data{j},2) > 1
                tempData.data(:,j) = ModelInfo.data{j}(:,i);
            else
                tempData.data(:,j) = ModelInfo.data{j};
            end
        end
        % Shuffle the specified column
        % tempData.data(:,tempData.ColumnToShuffle) = tempData.data(Samp(:,k),tempData.ColumnToShuffle);

        Results = FitProcessModelPermute(tempData,Samp(:,k));
        
        if Nperm == 0
            % this is the point estimate, keep everything
            Parameters{i} = Results;
        else
            % only keep the path estimates
            PermResults.Paths{i} = Results.Paths;
            PermResults.PathsTnorm{i} = Results.PathsTnorm;
            PermResults.beta(:,:,i) = Results.beta;
            PermResults.t(:,:,i) = Results.t;
            
        end
        
        % Need to calculate the cluster enhanced point estimate values!
        
        
        
        % Doing this is good to see progress but it SLOWS it down A LOT!!!
        % if ~mod(i,1000)
        %     fprintf(1,'Finished voxel %d of %d in %0.2f s.\n',i,Nvox,toc);
        % end
    end
    if Nperm > 0
        % fprintf(1,'\tFinding the maxima and minima from this permutaion.\n');
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
                    
                    % Need to add the TFCE calculations here also
                    % Perform all TFCE steps
                    % Take each t map and save it as a file
                    t = squeeze(ThisPath(j,m,:));
                    
                    % Z norm the path values so they are in the correct
                    % range for the TFCE algorithm
                    tempI(tempData.Indices) = t;
                    spm_write_vol(Vo,tempI);
%%%% THIS NEEDS TO BE RESHAPED %%%%%%%%%%%%%%%%5                    
                    [maxTFCE, minTFCE] = subfnApplyTFCE(Vo.fname, Vm.fname,TFCEparams);
                    
                    TFCEpathsMax(i,j,k) = maxTFCE;
                    TFCEpathsMin(i,j,k) = minTFCE;

                end
            end
        end
        % find the max and min beta values also
        for i = 1:Mbeta
            for j = 1:Nbeta
                if PermResults.beta(i,j,1) ~= 0
                    %fprintf(1,'%d, %d\n',i,j);
                    
                    MaxBeta(i,j,k) = max(squeeze(PermResults.beta(i,j,:)));
                    MinBeta(i,j,k) = min(squeeze(PermResults.beta(i,j,:)));
                    
                    % Perform all TFCE steps
                    % Take each t map and save it as a file
                    t = squeeze(PermResults.t(i,j,:,1));
                    tempI(tempData.Indices) = t;
                    spm_write_vol(Vo,tempI);
                    
                    [maxTFCE, minTFCE] = subfnApplyTFCE(Vo.fname, Vm.fname,TFCEparams);
                    
                    TFCEtMax(i,j,k) = maxTFCE;
                    TFCEtMin(i,j,k) = minTFCE;
                    
                end
            end
        end
        
         OutFile = fullfile(ResultsFolder,sprintf('Path_count%04d',k));
    Str = sprintf('save %s PermResults',OutFile);
    eval(Str)   
    end

    fprintf(1,'Finished permutation %d of %d in %0.2f s.\n',k,Nperm,toc);
end
fprintf(1,'Saving data to file now.\n\n');
% Save results to file
if Nperm > 0
    [PathName FileName] = fileparts(InDataPath);
    [PathName FileName] = fileparts(PathName);
    ResultsFolder = fullfile(PathName,'Results');
    if ~exist(ResultsFolder,'dir')
        mkdir(ResultsFolder)
    end
    OutFile = fullfile(ResultsFolder,sprintf('Permute_count%04d_%dSamp',count,Nperm));
    Str = sprintf('save %s MaxPaths MinPaths MaxBeta MinBeta TFCEtMax TFCEtMin TFCEpathsMax TFCEpathsMin',OutFile);
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
    fprintf(1,'Saved point estimate results to file.\n\n');
end
    

  
    
    

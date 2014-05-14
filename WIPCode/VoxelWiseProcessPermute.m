function VoxelWiseProcessPermute(InDataPath,count,Nperm)
% Add a check to determine if a path or a data structure is passed.

RandStream.setDefaultStream(RandStream('mt19937ar','Seed',sum(100*clock)));
fprintf(1,'Started at: %s\n',datestr(now));
tic
% load data
if nargin == 3
    % Load the 
    load(InDataPath)
    count = str2num(count);
    Nperm = str2num(Nperm);
elseif nargin == 2
    load(InDataPath)
    count = str2num(count);
    Nperm = 0;
else
    error('Expected two inputs');
end
NSub = ModelInfo.NSub;
if Nperm == 0
    % This is the point estimate
    Samp = [1:NSub]';
else
    Samp = zeros(NSub,Nperm);
    for i = 1:Nperm
        Samp(:,i) = randperm(NSub)';
    end
end
Nvox = length(ModelInfo.Indices);
%%
if Nperm > 0
    
    % create the structure to hold the permutation resample results
    % test one voxel to make sure the output structure is the correct size
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
    
    TestResults = WIPsubfnFitModel(tempData);
    % create the structure to hold the permutation results
    PermResults = {};
    FieldNames = {'Paths'};
    for k = 1:length(FieldNames)
        Value = getfield(TestResults,FieldNames{k});
        if iscell(Value)
            BlankValue = cell(size(Value,1),Nvox);
        else
            BlankValue = zeros([size(Value) Nvox]);
        end
        PermResults = setfield(PermResults,FieldNames{k},BlankValue);
    end
    
    [Mbeta Nbeta] = size(TestResults.B);
    PermResults.beta = zeros(Mbeta, Nbeta, Nvox, Nperm);
    
    [PathSize1 PathSize2] = size(TestResults.Paths{1});
    NPaths = length(TestResults.Paths);
    MaxB = zeros([size(TestResults.B) Nperm]);
    MinB = zeros([size(TestResults.B) Nperm]);
    MaxPaths = zeros(PathSize1,PathSize2,NPaths,Nperm);
    MinPaths = zeros(PathSize1,PathSize2,NPaths,Nperm);
end
%%


fprintf(1,'Data prepared in %0.2f s.\n',toc);
fprintf(1,'Starting first permutation at: %s\n',datestr(now));


for k = 1:size(Samp,2)
    tic
    for i = 1:Nvox
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
            ThisPath = zeros(PathSize1,PathSize2,Nvox);
            for j = 1:Nvox
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
    
    
    
    
    
    
    
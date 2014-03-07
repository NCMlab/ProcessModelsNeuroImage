function Parameters = WIPsubfnVoxelWiseProcessBatch(InData,BootStrapFlag)
% This program takes the data or the pointer to the data file and prepares
% the data for processing. It is here that multiple voxels are looped over.
% Therfore, it is in this program that one voxel at a time is extracted to
% create a temporary data array which is preocessed and the
% Parameters/results are returned. Therefore, if the voxelwise data is not
% the mediator then this program decides where the voxelwise data is and
% pulls out one voxel.
%
if nargin == 1
    BootStrapFlag = 0;
end
warning off

addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode

% Check to see if a path to data is given. If so then load it up
if ischar(InData)
    load(InData);
    % find out the tag for this batch of data
    [PathName FileName] = fileparts(InData);
    [PathName2 FileName2] = fileparts(PathName)
    ResultsFolder = fullfile(PathName2,'Results');
    if ~exist(ResultsFolder)
        mkdir(ResultsFolder);
    end
    
    tag = InData(end-3:end);
% if a structure is passed then just operate on this one data structure
elseif isstruct(InData)
    data = InData;
    Nboot = data.Nboot;
end
% Check the dimensions of the data
% The number of mediators can only be the second dimension of the mediator
% variable.
Nmed = size(data.M,2);
% The number of subjects will always be the first size dimension of any of
% the variables.
NSub = size(data.Y,1);
% take the number of voxels from the length of the indices instead of from
% the dimensions of the data because any of the variables could be
% voxel-wise.
Nvoxels = length(data.Indices);
% Initialize the cell array to hold the results
Parameters = cell(Nvoxels,1);

% Initialize the probe flag
data.ProbeMod = 0;

% Check to see if the variables are named
% If not use default names
if ~isfield(data.names,'X'); data.names.X = 'X'; end
if ~isfield(data.names,'M'); data.names.M = 'M'; end
if ~isfield(data.names,'Y'); data.names.Y = 'Y'; end
if ~isfield(data.names,'V'); data.names.V = 'V'; end
if ~isfield(data.names,'W'); data.names.W = 'W'; end
temp = data;
if ~isfield(data.names,'COV') && ~isempty(data.COV)
    NCov = size(data.COV,2);
    NameCovStruct = cell(NCov,1);
        for j = 1:NCov 
            NameCovStruct{j} = sprintf('Cov%d',j);
        end
    data.names.COV = NameCovStruct;    
    temp.names.COV = NameCovStruct;
end
Samp = 1:NSub;
tic 
for i = 1:Nvoxels

    % check to make sure there is data for all subjects at this voxel. 
    Mflag = 0;
    Vflag = 0;
    Wflag = 0;
% added this flag so that the model check cases do the checking and then
% they set this flag
    AllDataFlag = 0;
    
    switch data.ModelNum
        case '1'
            % Set the probeMode flag to TRUE for the first call 
            if Nmed ~= 1
                errordlg('Only a single variable can be used as a moderator for Model 1');
            end
            temp = CheckVariables(i,data,temp);
           
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)&(sum(isnan(temp.Y)) == 0)
                AllDataFlag = 1;
            end

        case '4'
            % Check each variable and if it is a voxel-wise measure then
            % pull out one voxel. It is also possible to have all varaibles
            % be voxelwise and to have voxelwise covariates.
            % Check the X variable
            temp = CheckVariables(i,data,temp);
            % If there are no not-a-number variables then the data is ready
            % to be processed.
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)&(sum(isnan(temp.Y)) == 0)
                AllDataFlag = 1;
            end
            
        case '6'
           
            temp = CheckVariables(i,data,temp);
            tempData = {};
            tempData.data = [temp.X temp.M temp.Y];
            tempData.STRAT = temp.STRAT;
            M = size(tempData.data,2);
            Direct = zeros(M,M);
            Direct(1,[2 3 4]) = 1;
            Direct(2,[3 4]) = 1;
            Direct(3,[4]) = 1;
            Inter = zeros(M,M);
            %Inter([1 3],4) = 1;
            Paths = zeros(M,M,3);
            Paths(1,2,1) = 1;
            Paths(2,3,1) = 2;
            Paths(3,4,1) = 3;
            Paths(1,2,2) = 1;
            Paths(2,4,2) = 2;
            Paths(1,3,3) = 1;
            Paths(3,4,3) = 2;
            tempData.Direct = Direct;
            tempData.Inter = Inter;
            tempData.Paths = Paths;

            % If there are no not-a-number variables then the data is ready
            % to be processed.
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)&(sum(isnan(temp.Y)) == 0)
                AllDataFlag = 1;
            end
        case '7'
            if isempty(data.W)
                errordlg('The modulator variable is missing');
            end
            % Calculate the minimum critical t-value for use with the
            % Johnson-Neyman technique.
            temp.tcrit = tinv(1 - max(data.Thresholds),NSub - 4);
            % Check each variable and if it is a voxel-wise measure then
            % pull out one voxel. It is also possible to have all varaibles
            % be voxelwise and to have voxelwise covariates.
            % Check the X variable
            temp = CheckVariables(i,data,temp);
            % If there are no not-a-number variables then the data is ready
            % to be processed.
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)...
                    &(sum(isnan(temp.Y)) == 0)&(sum(isnan(temp.W)) == 0)
                AllDataFlag = 1;
            end
        case '14'
            if isempty(data.V)
                errordlg('The modulator variable is missing');
            end
            
            % Calculate the minimum critical t-value for use with the
            % Johnson-Neyman technique.
            temp.tcrit = tinv(1 - max(data.Thresholds),NSub - 4);
            temp = CheckVariables(i,data,temp);
            % If there are no not-a-number variables then the data is ready
            % to be processed.
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)...
                    &(sum(isnan(temp.Y)) == 0)&(sum(isnan(temp.V)) == 0)
                AllDataFlag = 1;
            end

        case '58'
            if isempty(data.W)
                errordlg('The modulator variable is missing');
            end
            % Calculate the minimum critical t-value for use with the
            % Johnson-Neyman technique.
            temp.tcrit = tinv(1 - max(data.Thresholds),NSub - 4);
            temp = CheckVariables(i,data,temp);
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)...
                    &(sum(isnan(temp.Y)) == 0)...
                    &(sum(isnan(temp.W)) == 0)
                AllDataFlag = 1;
            end
          
        case '74'
            temp.tcrit = tinv(1 - max(data.Thresholds),NSub - 4);
            temp = CheckVariables(i,data,temp);
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)&(sum(isnan(temp.Y)) == 0)
                AllDataFlag = 1;
            end
            
        case '75'
            temp = CheckVariables(i,data,temp);
            tempData = {};
            tempData.data = [temp.X temp.M temp.Y];
            tempData.STRAT = temp.STRAT;
            M = size(tempData.data,2);
            Direct = zeros(M,M);
            Direct(1,[2 3 4]) = 1;
            Direct(2,[3 4]) = 1;
            Direct(3,[4]) = 1;
            Inter = zeros(M,M);
            Inter([1 3],4) = 1;
            Paths = zeros(M,M,1);
            Paths(1,2) = 1;
            Paths(2,3) = 2;
            Paths(3,4) = 3;
            tempData.Direct = Direct;
            tempData.Inter = Inter;
            tempData.Paths = Paths;
            % If there are no not-a-number variables then the data is ready
            % to be processed.
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)&(sum(isnan(temp.Y)) == 0)
                AllDataFlag = 1;
            end
    end
    % Check for data being all zeros
    if length(find(temp.X==0))==length(temp.X)
        AllDataFlag = 0;
    end
    if length(find(temp.Y==0))==length(temp.Y)
        AllDataFlag = 0;
    end
    if length(find(temp.M==0))==length(temp.M)
        AllDataFlag = 0;
    end
    if AllDataFlag
        if BootStrapFlag
            % create the resamples only once
            Samp = BootStrap(tempData);
            BootStrapFlag = 0;
        end
        tempData.data = tempData.data(Samp,:);
        
        Results = WIPsubfnFitModel(tempData);
        Parameters{i} = Results;
        %Parameters{i}.Nboot = data.Nboot;
        %Parameters{i}.Thresholds = data.Thresholds;
        %fprintf(1,'\n')
    else
        fprintf(1,' << empty\n')
    end
%     if Nvoxels > 1
%         if ~mod(i,1000)
%             t = toc;
%             fprintf(1,'%6d of %6d in %6.2f sec \n',i,Nvoxels,t);
%             tic;
%         end
%     end
end

if ischar(InData)
    % save the results to a mat file so that the main program can load them up
    Str = ['save ' fullfile(ResultsFolder,['Results_' tag]) ' Parameters'];
    eval(Str);
    fprintf(1,'Done!');
end


function Samp = BootStrap(data)
% For each cluster chunk it is creating the EXACT same set of bootstrap
% resamples ... this is bad. The following line attempts to avoid that by
% resetting the random seed.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
if ~isempty(data.STRAT)
    Gr1 = find(data.STRAT == 0);
    NGr1 = length(Gr1);
    Gr2 = find(data.STRAT == 1);
    NGr2 = length(Gr2);
else
    NGr1 = [];
    NGr2 = [];
end
temp = data;
% create the resamples
if isempty(data.STRAT)
    N = size(data.data,1);
    Samp =  floor(N*rand(N,1))+1;
else
    Samp1 = floor(NGr1*rand(NGr1,1))+1;
    Samp2 = floor(NGr2*rand(NGr2,1))+1+NGr1;
    Samp = [Samp1; Samp2];
end




function temp = CheckVariables(VoxelIndex,data,temp)
if size(data.X,3) > 1
    temp.X = data.X(:,:,VoxelIndex);
else
    temp.X = data.X;
end
% Check the M variable
if size(data.M,3) > 1
    temp.M = data.M(:,:,VoxelIndex);
else
    temp.M = data.M;
end
% Check the Y variable
if size(data.Y,3) > 1
    temp.Y = data.Y(:,:,VoxelIndex);
else
    temp.Y = data.Y;
end
% Check the V variable
if size(data.V,3) > 1
    temp.V = data.V(:,:,VoxelIndex);
else
    temp.V = data.V;
end
% Check the W variable
if size(data.W,3) > 1
    temp.W = data.W(:,:,VoxelIndex);
else
    temp.W = data.W;
end
% Check the Q variable
if size(data.Q,3) > 1
    temp.Q = data.Q(:,:,VoxelIndex);
else
    temp.Q = data.Q;
end
% Check the R variable
if size(data.R,3) > 1
    temp.R = data.R(:,:,VoxelIndex);
else
    temp.R = data.R;
end
% Check the covariates
if size(data.COV,3) > 1
    temp.COV = data.COV(:,VoxelIndex);
else
    temp.COV = data.COV;
end

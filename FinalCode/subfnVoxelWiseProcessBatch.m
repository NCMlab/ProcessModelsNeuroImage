function Parameters = subfnVoxelWiseProcessBatch(InData)
% This program takes the data or the pointer to the data file and prepares
% the data for processing. It is here that multiple voxels are looped over.
% Therfore, it is in this program that one voxel at a time is extracted to
% create a temporary data array which is preocessed and the
% Parameters/results are returned. Therefore, if the voxelwise data is not
% the mediator then this program decides where the voxelwise data is and
% pulls out one voxel.
%
warning off
addpath /share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode
%
% Check to see if a path to data is given. If so then load it up
if ischar(InData)
    load(InData);
    % find out the tag for this batch of data
    [PathName FileName] = fileparts(InData);
    tag = InData(end-3:end);
    [Nsub] = size(data.M,1);
    Nmed = size(data.M,2);
    Nvoxels = length(data.Indices);
    Parameters = cell(Nvoxels,1);

    % If this voxelwise data than try to run in parallel
    
%     try
%         NumOpen = matlabpool('size');
%         MaxPool = 8;
%         if NumOpen < MaxPool
%             Str = ['matlabpool open ' num2str(MaxPool - NumOpen)];
%             eval(Str)
%         end
%         OpenPoolFlag = 1;
%     catch
%         OpenPoolFlag = 0;
%     end
    
    
% if a structure is passed then just operate on this one data structure
elseif isstruct(InData)
    data = InData;
    [Nsub Nmed Nvoxels] = size(data.M);
    Parameters = cell(Nvoxels,1);
    Nboot = data.Nboot;
end
data.ProbeMod = 0;

% Check to see if the variables are named
% If not use default names
if ~isfield(data,'Xname'); data.Xname = 'X'; end
if ~isfield(data,'Mname'); data.Mname = 'M'; end
if ~isfield(data,'Yname'); data.Yname = 'Y'; end
if ~isfield(data,'Vname'); data.Vname = 'V'; end
if ~isfield(data,'Wname'); data.Wname = 'W'; end
temp = data;
if ~isfield(data,'COVname') && ~isempty(data.COV)
    NCov = size(data.COV,2);
    NameCovStruct = cell(NCov,1);
        for j = 1:NCov 
            NameCovStruct{j} = sprintf('Cov%d',j);
        end
    data.COVname = NameCovStruct;    
    temp.COVname = NameCovStruct;
end

tic 
for i = 1:Nvoxels
    if Nvoxels > 1
        t = toc;
        fprintf(1,'%6d of %6d in %6.2f sec\n',i,Nvoxels,t);
        tic;
    end
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
            if size(data.M,2) ~= 1
                errordlg('Only a single variable can be used as a moderator for Model 1');
            end
            temp.M = data.M(:,:,i);
            if sum(isnan(temp.M)) == 0
                AllDataFlag = 1;
            end

        case '4'
            % Check each variable and if it is a voxel-wise measure then
            % pull out one voxel. It is also possible to have all varaibles
            % be voxelwise and to have voxelwise covariates.
            % Check the X variable
            if size(data.X,3) > 1
                temp.X = data.X(:,:,i);
            else
                temp.X = data.X;
            end
            % Check the M variable
            if size(data.M,3) > 1
                temp.M = data.M(:,:,i);
            else
                temp.M = data.M;
            end
            % Check the Y variable
            if size(data.Y,3) > 1
                temp.Y = data.Y(:,:,i);
            else
                temp.Y = data.Y;
            end
            % Check the covariates
            if size(data.COV,3) > 1
                temp.COV = data.COV(:,i);
            else
                temp.COV = data.COV;
            end
            % If there are no not-a-number variables then the data is ready
            % to be processed.
            if (sum(isnan(temp.X)) == 0)&(sum(isnan(temp.M)) == 0)&(sum(isnan(temp.Y)) == 0)
                AllDataFlag = 1;
            end
        case '7'
            if isempty(data.W)
                errordlg('The modulator variable is missing');
            end
            if size(data.M,2) ~= 1
                errordlg('Only a single variable can be used as a mediator for Model 7');
            end
            temp.M = data.M(:,:,i);
            temp.W = data.W(:,:,i);
            if sum(isnan(temp.M)) == 0; Mflag = 1;end
            if sum(isnan(temp.W)) == 0; Wflag = 1;end

            if Mflag && Wflag 
                AllDataFlag = 1;
            end
        case '14'
            if isempty(data.V)
                errordlg('The modulator variable is missing');
            end
            % Calculate the minimum critical t-value for use with the
            % Johnson-Neyman technique.
            temp.tcrit = tinv(1 - max(data.Thresholds),Nsub - 4);
            temp.M = data.M(:,:,i);
            temp.V = data.V(:,:,i);
            if sum(isnan(temp.M)) == 0; Mflag = 1;end
            if sum(isnan(temp.V)) == 0; Vflag = 1;end

            if Mflag && Vflag 

                AllDataFlag = 1;
            end
    end
    
    if AllDataFlag
        tempParameters = subfnProcess(temp);
        Parameters{i} = tempParameters{:};
        Parameters{i}.Nboot = data.Nboot;
        Parameters{i}.Thresholds = data.Thresholds;
    end

end

if ischar(InData)
    % save the results to a mat file so that the main program can load them up
    Str = ['save ' fullfile(PathName,['Results_' tag]) ' Parameters'];
    eval(Str);
    fprintf(1,'Done!')
end

% if OpenPoolFlag
%     matlabpool close 
% end


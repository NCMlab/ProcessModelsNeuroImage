function Parameters = subfnVoxelWiseProcessBatch(InData)
addpath /share/data/users/js2746_Jason/CogReserveAnalyses/FinalCode
% OpenPoolFlag = 1;
% try
%     NumOpen = matlabpool('size');
%     MaxPool = 8;
%     if NumOpen < MaxPool
%         Str = ['matlabpool open ' num2str(MaxPool - NumOpen)];
%         eval(Str)
%     end
% catch 
%     OpenPoolFlag = 0;
% end
%
% Check to see if a path to data is given. If so then load it up
if ischar(InData)
    load(InData);
    % find out the tag for this batch of data
    [PathName FileName] = fileparts(InData);
    tag = InData(end-3:end);
    [Nsub Nmed Nvoxels] = size(data.M);
    temp = data;
    Parameters = cell(Nvoxels,1);
% if a structure is passed then just operate on this one data structure
elseif isstruct(InData)
    temp = InData;
    data = temp;
    [Nsub Nmed Nvoxels] = size(temp.M);
    Parameters = cell(Nvoxels,1);
    ModelNum = temp.ModelNum;
    Nboot = temp.Nboot;
    Thresholds = temp.Thresholds;
end

for i = 1:Nvoxels
    % check to make sure there is data for all subjects at this voxel
    Mflag = 0;
    Vflag = 0;
    switch ModelNum
        case '4'
            M = data.M(:,:,i);
            V = [];
            if sum(isnan(M)) == 0
                Mflag = 1;
                Vflag = 1;
            end
        case '14'
            M = data.M(:,:,i);
            V = data.V(:,:,i);
            if sum(isnan(M)) == 0; Mflag = 1;end
            if sum(isnan(V)) == 0; Vflag = 1;end
    end
    if Mflag & Vflag
        temp.M = M;
        temp.V = V;
        
        [pointEst Parameters{i}] = subfnProcessModelFit(temp,ModelNum,1);
        % have the option of turing the boot strapping off
        if Nboot
            % calculate the boot strap estimates
            [bstat] = subfnBootStrp(temp,Nboot,ModelNum);
            % calculate the BCAlimits
            [BCaci PERci] = subfnFindConfidenceIntervals(temp,bstat,pointEst,ModelNum,Thresholds);
            for j = 1:size(bstat,2)
                for k = 1:size(bstat,3)
                    % there is an issue here with the
                    %Parameters{i}.AB{k,j}.pointEst = pointEst(k);
                    Parameters{i}.AB{k,j}.BCaci = BCaci{j,k};
                    Parameters{i}.AB{k,j}.PERci = PERci{j,k};
                    Parameters{i}.AB{k,j}.se = bootSE(bstat(:,j),Nsub);
                end
                Parameters{i}.Nboot = Nboot;
                Parameters{i}.VoxelIndex = data.Indices(i);
            end
        else
            for j = 1:size(pointEst,1);
                for k = 1:size(pointEst,2);
                    Parameters{i}.AB{k,j}.BCaci = [];
                    Parameters{i}.AB{k,j}.PERci = [];
                    Parameters{i}.AB{k,j}.se = [];
                end
            end
        end
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
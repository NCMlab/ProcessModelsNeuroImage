function Parameters = subfnVoxelWiseProcessBatch(InData)
addpath /share/data/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode

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
    % check to make sure there is data for all subjects at this voxel. 
    Mflag = 0;
    Vflag = 0;
    temp.ProbeMod = 0;
    switch ModelNum
        case '1'
            % Set the probeMode flag to TRUE for the first call 
           
            M = data.M(:,:,i);
            V = [];
            if sum(isnan(M)) == 0
                Mflag = 1;
                Vflag = 1;
            end
        case '4'
            M = data.M(:,:,i);
            V = [];
            if sum(isnan(M)) == 0
                Mflag = 1;
                Vflag = 1;
            end
        case '14'
            % Calculate the minimum critical t-value for use with the
            % Johnson-Neyman technique.
            temp.tcrit = tinv(1 - max(Thresholds),Nsub - 4);
            M = data.M(:,:,i);
            V = data.V(:,:,i);
            if sum(isnan(M)) == 0; Mflag = 1;end
            if sum(isnan(V)) == 0; Vflag = 1;end
    end
    
    if Mflag & Vflag
        temp.M = M;
        temp.V = V;
        % it initially does not probe the moderator
        [pointEst Parameters{i}] = subfnProcessModelFit(temp,1);
        % here we check to see if the interaction is significant.
        % if so then probe the mod for all bootstraps
        if isfield(Parameters{i},'Int')
            if Parameters{i}.Int{1}.p < max(data.Thresholds)
                temp.ProbeMod = 1;
                % since the interaction is significant then we want to probe
                % the moderator by re-running the regression
               [pointEst Parameters{i}] = subfnProcessModelFit(temp,1);
            end
        end
        % have the option of turning the boot strapping off.
        % turn off boot strapping on the interaction effect if the
        % Johnson-Neyman shows no range of significant interactions.
        %
        if Nboot %& Parameters{i}.JohnsonNeyman ~= -99 
            % calculate the boot strap estimates
            [bstat] = subfnBootStrp(temp,Nboot);
            [nboot Nmed NParameters] = size(bstat);
            % Find the number of non-zero bstat values
            if length(find(squeeze(bstat(1,1,:)))) == 1
                NParameters = 1;
            end
            
            % calculate the Confidence intervals on the parameters that
            % need to be bootstrapped.
            [BCaci PERci] = subfnFindConfidenceIntervals(temp,bstat,pointEst,Thresholds);
            % Then fill in the appropriate Parameters with the confidence
            % intervals.
            str = [pointEst.names ' = {};'];
            eval(str)

            for j = 1:Nmed
                temp2 = {};
                for k = 1:NParameters
%                     str = sprintf('temp2 = Parameters{i}.%s{j};',pointEst.names);
%                     eval(str);
                    str = sprintf('temp2=setfield(temp2,''BCaci'',BCaci{%d,%d});',j,k);
                    eval(str);
                    str = sprintf('temp2=setfield(temp2,''PERci'',PERci{%d,%d});',j,k);
                    eval(str);
                    str = sprintf('temp2=setfield(temp2,''bootSE'',bootSE(bstat(:,%d,%d),Nsub));',j,k);
                    eval(str);
                    str = sprintf('temp2=setfield(temp2,''pointEst'',pointEst.values(k));');
                    eval(str);
                    if isfield(pointEst,'probeValues')
                        str = sprintf('temp2=setfield(temp2,''probeValue'',pointEst.probeValues(k));');
                        eval(str);
                    end
                    str = sprintf('%s{%d,%d}=temp2;',pointEst.names,k,j);
                    eval(str);
                end
                str = sprintf('Parameters{i}.%s{j} = %s;',pointEst.names,pointEst.names);
                eval(str);
            end
            
        else
            for j = 1:size(pointEst,1);
                for k = 1:size(pointEst,2);
                    str = sprintf('temp2=setfield(temp2,''BCaci'',[]);');
                    eval(str);
                    str = sprintf('temp2=setfield(temp2,''PERci'',[]);');
                    eval(str);
                    str = sprintf('temp2=setfield(temp2,''bootSE'',[]);');
                    eval(str);
                    str = sprintf('%s{%d,%d}=temp2;',pointEst.names,k,j);
                    eval(str);
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
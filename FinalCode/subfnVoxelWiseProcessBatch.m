function Parameters = subfnVoxelWiseProcessBatch(InData)
<<<<<<< HEAD
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
=======

warning off
addpath /share/data/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode

OpenPoolFlag = 0;
>>>>>>> develop
%
% Check to see if a path to data is given. If so then load it up
if ischar(InData)
    load(InData);
    % find out the tag for this batch of data
    [PathName FileName] = fileparts(InData);
    tag = InData(end-3:end);
    [Nsub Nmed Nvoxels] = size(data.M);
<<<<<<< HEAD
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
    Wflag = 0;
    % added this flag so that the model check cases do the checking and then
    % they set this flag
    AllDataFlag = 0;
    temp.ProbeMod = 0;
    switch ModelNum
        case '1'
            % Set the probeMode flag to TRUE for the first call
=======
    Parameters = cell(Nvoxels,1);
    data.Thresholds = Thresholds;
    data.ModelNum = ModelNum;
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
    t = toc;
    fprintf(1,'%6d of %6d in %6.2f sec\n',i,Nvoxels,t);
    tic;
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
>>>>>>> develop
            if size(data.M,2) ~= 1
                errordlg('Only a single variable can be used as a moderator for Model 1');
            end
            temp.M = data.M(:,:,i);
            if sum(isnan(temp.M)) == 0
                AllDataFlag = 1;
            end
<<<<<<< HEAD
            
        case '4'
            
            temp.M = data.M(:,:,i);
            if sum(isnan(temp.M)) == 0;
=======

        case '4'
            
            temp.M = data.M(:,:,i);
            if sum(isnan(temp.M)) == 0; 
>>>>>>> develop
                AllDataFlag = 1;
            end
            
        case '7'
            if isempty(data.W)
                errordlg('The modulator variable is missing');
            end
            temp.M = data.M(:,:,i);
            temp.W = data.W(:,:,i);
            if sum(isnan(temp.M)) == 0; Mflag = 1;end
            if sum(isnan(temp.W)) == 0; Wflag = 1;end
<<<<<<< HEAD
            if Mflag && Wflag
=======
            if Mflag && Wflag 
>>>>>>> develop
                AllDataFlag = 1;
            end
        case '14'
            if isempty(data.V)
                errordlg('The modulator variable is missing');
            end
            % Calculate the minimum critical t-value for use with the
            % Johnson-Neyman technique.
            temp.tcrit = tinv(1 - max(Thresholds),Nsub - 4);
            temp.M = data.M(:,:,i);
            temp.V = data.V(:,:,i);
            if sum(isnan(temp.M)) == 0; Mflag = 1;end
            if sum(isnan(temp.V)) == 0; Vflag = 1;end
<<<<<<< HEAD
            if Mflag && Vflag
=======
            if Mflag && Vflag 
>>>>>>> develop
                AllDataFlag = 1;
            end
    end
    
    if AllDataFlag
        % it initially does not probe the moderator
        [pointEst Parameters{i}] = subfnProcessModelFit(temp,1);
        % here we check to see if the interaction is significant.
        % if so then probe the mod for all bootstraps
        if pointEst.probeMod
<<<<<<< HEAD
            temp.ProbeMod = 1;
            % since the interaction is significant then we want to probe
            % the moderator by re-running the regression
            [pointEst Parameters{i}] = subfnProcessModelFit(temp,1);
=======
                temp.ProbeMod = 1;
                % since the interaction is significant then we want to probe
                % the moderator by re-running the regression
               [pointEst Parameters{i}] = subfnProcessModelFit(temp,1);
>>>>>>> develop
            
        end
        % have the option of turning the boot strapping off.
        % turn off boot strapping on the interaction effect if the
        % Johnson-Neyman shows no range of significant interactions.
        %
<<<<<<< HEAD
        if Nboot %& Parameters{i}.JohnsonNeyman ~= -99
            % calculate the boot strap estimates
            [bstat] = subfnBootStrp(temp,Nboot);
=======
        if Nboot %& Parameters{i}.JohnsonNeyman ~= -99 
            % calculate the boot strap estimates
            [bstat k2stat] = subfnBootStrp(temp,Nboot);
            if k2stat(1) ~= 0
                Sk2stat = sort(k2stat);
                k2 = {};
                k2.pointEst = pointEst.k2;
                PERci = cell(length(temp.Thresholds));
                
                for k = 1:length(temp.Thresholds)
                    t = num2str(temp.Thresholds(k));
                    PERci{k} = setfield(PERci{k},['alpha' t(3:end)],...
                        [Sk2stat(ceil(temp.Thresholds(1)/2*Nboot)) Sk2stat(ceil((1-temp.Thresholds(1)/2)*Nboot))]);
                end
                k2.PERci = PERci;
                Parameters{i}.k2 = k2;
            end
            
>>>>>>> develop
            [nboot Nmed NParameters] = size(bstat);
            % Find the number of non-zero bstat values
            if length(find(squeeze(bstat(1,1,:)))) == 1
                NParameters = 1;
            end
            
            % calculate the Confidence intervals on the parameters that
            % need to be bootstrapped.
<<<<<<< HEAD
            [BCaci PERci] = subfnFindConfidenceIntervals(temp,bstat,pointEst,Thresholds);
            % Then fill in the appropriate Parameters with the confidence
            % intervals.
            %            str = [pointEst.names ' = {};'];
            %            eval(str)
            
            for j = 1:Nmed
                temp2 = {};
                for k = 1:NParameters
                    %                     str = sprintf('temp2 = Parameters{i}.%s{j};',pointEst.names);
                    %                     eval(str);
=======
            [BCaci PERci] = subfnFindConfidenceIntervals(temp,bstat,pointEst);
            % Then fill in the appropriate Parameters with the confidence
            % intervals.
%            str = [pointEst.names ' = {};'];
%            eval(str)

            for j = 1:Nmed
                temp2 = {};
                for k = 1:NParameters
%                     str = sprintf('temp2 = Parameters{i}.%s{j};',pointEst.names);
%                     eval(str);
>>>>>>> develop
                    str = sprintf('temp2=setfield(temp2,''BCaci'',BCaci{%d,%d});',j,k);
                    eval(str);
                    str = sprintf('temp2=setfield(temp2,''PERci'',PERci{%d,%d});',j,k);
                    eval(str);
                    str = sprintf('temp2=setfield(temp2,''bootSE'',bootSE(bstat(:,%d,%d),Nsub));',j,k);
                    eval(str);
                    str = sprintf('temp2=setfield(temp2,''pointEst'',pointEst.values(j,k));');
                    eval(str);
                    if isfield(pointEst,'probeValues')
                        str = sprintf('temp2=setfield(temp2,''probeValue'',pointEst.probeValues(k));');
                        eval(str);
                    end
                    str = sprintf('%s{%d}=temp2;',pointEst.names(j,:),k);
                    eval(str);
                end
                
                str = sprintf('Parameters{i}.%s = %s;',pointEst.names(j,:),pointEst.names(j,:));
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
<<<<<<< HEAD
% if OpenPoolFlag
%     matlabpool close
% end
=======
if OpenPoolFlag
    matlabpool close 
end
>>>>>>> develop

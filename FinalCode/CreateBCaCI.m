function BCaCI = CreateBCaCI(Results,BootStrap,JackKnife,Thresholds)
% from the point estimate, the bootstrap and the jacknife results return
% all of the bias-corrected, accelerated confidence intervals.

% Therefore, the whatever bootstrap estimates are calculated this function
% creates the BCa confidence intervals for it.

% Find the parameters in the bootstrap
FieldNames = fieldnames(BootStrap);

Nthresh = length(Thresholds);
NPaths = size(Results.Paths,1);
% prepare the structure to hold all of the BCa confidence intervals
BCaCI = {};
for i = 1:length(FieldNames)
    Value = getfield(Results,FieldNames{i});
    if iscell(Value)
        BlankValue = zeros(numel(Results.Paths{1}),1,2,NPaths,Nthresh);
    else
        BlankValue = zeros([size(Value) 2 Nthresh]);
    end
    BCaCI = setfield(BCaCI,FieldNames{i},BlankValue);
end
%% put the confidence intervals into this structure
% cycle over the input fieldnames
for i = 1:length(FieldNames)
    if ~isempty(strmatch(FieldNames{i},'Paths'))
        CurrentBCaCI = zeros(size(getfield(BCaCI,FieldNames{i})));
        % cycle over all the input thresholds
        for t = 1:Nthresh
            CurrentAlpha = Thresholds(t);
            for j = 1:NPaths
                % This takes that date from each path and flattens it
                BootStrapData = ConvertPathsToMatrix(BootStrap.Paths,j);
                JackKnifeData = ConvertPathsToMatrix(JackKnife.Paths,j);
                
                PointEstimate = reshape(Results.Paths{j},numel(Results.Paths{j}),1);
                
                % The flattended data allows it to be used by this
                % program which calculates the adjusted alpha limits used for
                % determining the confidence intervals
                [Alpha1 Alpha2 Z p] = CalculateBCaLimits(JackKnifeData,PointEstimate, BootStrapData,CurrentAlpha);
                
                % find the confidence intervals for these adjusted
                % alpha limits
                CurrentBCaCI(:,:,:,j,t) = CalculateBCaCI(BootStrapData,Alpha1,Alpha2,PointEstimate);
                BCaCI.PathsZ(j,:) = Z;
                BCaCI.PathsP(j,:) = p;

            end
        end
        % put the confidence interval data back into the structure
        % TO DO: somehow UNFLATTEN the data
        BCaCI = setfield(BCaCI,FieldNames{i},CurrentBCaCI);
    else
        CurrentBCaCI = zeros(size(getfield(BCaCI,FieldNames{i})));
        for t = 1:Nthresh
            CurrentAlpha = Thresholds(t);
            BootStrapData = getfield(BootStrap,FieldNames{i});
            JackKnifeData = getfield(JackKnife,FieldNames{i});
            PointEstimate = getfield(Results,FieldNames{i});
            [Alpha1 Alpha2 Z p] = CalculateBCaLimits(JackKnifeData,PointEstimate, BootStrapData,CurrentAlpha);
          
            % find the confidence intervals
            CurrentBCaCI(:,:,:,t) = CalculateBCaCI(BootStrapData,Alpha1,Alpha2,PointEstimate);
        end
        BCaCI = setfield(BCaCI,FieldNames{i},CurrentBCaCI);
        BCaCI.Z = Z;
        BCaCI.p = p;
    end
end
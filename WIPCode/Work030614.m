

data = {};
data{1} = AgeGroup1old;
data{2} = squeeze(STRUCTURE);
data{3} = squeeze(FMRI);
data{4} = BEHAVIOR;
Names = {};
Names{1} = 'AgeGr';
Names{2} = 'pGM';
Names{3} = 'fMRI';
Names{4} = 'medRTSC';


Nvar = length(data);
Direct = zeros(Nvar);
Inter = zeros(Nvar);

Paths = zeros(Nvar,Nvar,1);

Direct(1,[2 3 4]) = 1;
Direct(2,[3 4]) = 1;
Direct(3,[4]) = 1;
Inter([1 3],4) = 1;
Paths(1,2,1) = 1;
Paths(2,3,1) = 2;
Paths(3,4,1) = 3;

Nboot = 5000;
Thresh = [0.05 0.01 0.001];

ModelInfo = {};
ModelInfo.Names = Names;
ModelInfo.data = data;
ModelInfo.Direct = Direct;
ModelInfo.Inter = Inter;
ModelInfo.Paths = Paths;
ModelInfo.Nboot = Nboot;
ModelInfo.BaseDir = BaseDir;
ModelInfo.Tag= 'JASON';
ModelInfo.Indices = Indices;
ModelInfo.NJobSplit = 100;
ModelInfo.Thresholds = Thresh;
ModelInfo.STRAT = [];
ModelInfo.NSub = size(data{1},1);
ModelInfo.Nvar = Nvar;

% The following is what I am preparing as the main computational piece of
% code that will be sent out to the cluster nodes.
Nvox = length(data.Indices);
tempData = data;
Parameters = cell(Nvox,1);

FieldNames = {'beta' 'B' 'Paths'}
for i = 1:Nvox
    tic
    % for this voxel
    tempData.data = zeros(data.NSub,data.Nvar);
    % extract the data 
    for j = 1:data.Nvar
        if size(data.data{j},2) > 1
            tempData.data(:,j) = data.data{j}(:,i);
        else
            tempData.data(:,j) = data.data{j};
        end
    end
    % calculate the point estimate
    PointEstimateResults = WIPsubfnFitModel(tempData);
    Parameters{i} = PointEstimateResults;
    
    % Bootstrap this voxel
    BootStrap = {};
    % create the structure to hol dthe bootstrap resample results
    for k = 1:length(FieldNames)
        Value = getfield(Results,FieldNames{k});
        if iscell(Value)
            BlankValue = cell(size(Value),Nboot);
        else
            BlankValue = zeros([size(Value) Nboot]);
        end
        BootStrap = setfield(BootStrap,FieldNames{k},BlankValue);
    end
    % run the boot strapping
    for j = 1:data.Nboot
        BStempData = tempData;
        BStempData.data = tempData.data(data.BootSamp(:,j),:);
        BSResults = WIPsubfnFitModel(BStempData);
        BootStrap.beta(:,:,j) = BSResults.beta;
        BootStrap.B(:,:,j) = BSResults.B;
        BootStrap.Paths{j} = BSResults.Paths;
    end
    % Calculate the statistics on the BS results.
    toc
    % I am thinking it might be a good idea to save the BCa parameters like
    % the theta value from the jack-knife step.
    % Jack-knife
    Value = getfield(PointEstimateResults,'Paths');
    JackKnifeValue = cell(size(Value),NSub);
    JKtempData = tempData;
    for k = 1:NSub
        Include = ones(NSub,1);
        Include(k) = 0;
        JKtempData.data = tempData.data(logical(Include),:);
        JKResults = WIPsubfnFitModel(BStempData);
        JackKnifeValue{k} = JKResults.Paths;
    end
 %%   
Alpha1 = zeros(Nmed,NParameters);
Alpha2 = zeros(Nmed,NParameters);
p = zeros(Nmed,NParameters);
Z = zeros(Nmed,NParameters);


% 
testValue = 0;
for k = 1:NParameters
    for j = 1:Nmed
        zh0 = norminv(length(find(bstat(:,j,k) < pointEst(j,k)))/nboot);
        zA = norminv(alpha/2);
        z1mA = norminv(1 - alpha/2);
        ThetaDiff = (sum(theta(:,j))/N) - theta(:,j);
        % Calculate the acceleration factor
        acc = (sum(ThetaDiff.^3))/(6*(sum(ThetaDiff.^2))^(3/2));
        Alpha1(j,k) = normcdf(zh0 + (zh0+zA)/(1 - acc*(zh0 + zA)));
        Alpha2(j,k) = normcdf(zh0 + (zh0+z1mA)/(1 - acc*(zh0 + z1mA)));
        % Find percentile of the distribution below the null value
        PCTlower = sum(bstat(:,j,k) < testValue)./nboot;
        PCTupper = sum(bstat(:,j,k) > testValue)./nboot;
        IsLower = PCTlower < PCTupper;
        
        
        % Check to make sure the tails are not zero
        PCTlower = max(PCTlower,1./nboot);
        PCTupper = min(PCTupper, 1-1./nboot);
        
        PCTupper = max(PCTupper,1/nboot);
        PCTlower = min(PCTlower,1-1/nboot);
        
        % Z-score
        ZPCTlower = norminv(PCTlower) - zh0;
        ZPCTupper = norminv(PCTupper) - zh0;
        % adjust for acceleration
        Zlower = (ZPCTlower.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTlower);
        Zupper = (ZPCTupper.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTupper);
        % which tail are we in?
        
        ZTemp = Zupper;
        ZTemp(IsLower) = -Zlower(IsLower);
        Z(j,k) = ZTemp;
        
        pTemp = normcdf(ZTemp);
        pTemp = [pTemp; 1 - pTemp];
        p(j,k) = 2.*min(pTemp);
        
        
    end
end
    
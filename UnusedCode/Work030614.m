addpath /share/studies/CogRes/Scripts/ProcessModelsNeuroImage/FinalCode
clear
BaseDir = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/SampleSize140';
cd(BaseDir)
load('Age_GM_ECF_PERF_LargerSample');

Nsub = length(AgeGroup1old);

Pmask = fullfile(BaseDir,'mask','maskIncGrey50.nii');
Vmask = spm_vol(Pmask);
Imask = spm_read_vols(Vmask);
Indices = find(Imask);
Nvoxels = length(Indices);

% Load up Functional Data
fprintf(1,'Loading fMRI data...\n');
FMRI = zeros(Nsub,1,Nvoxels);
for i = 1:Nsub
    V = spm_vol(FMRIPATHS{i});
    I = spm_read_vols(V);
    FMRI(i,1,:) = I(Indices);
end

% Load up Structural Data
fprintf(1,'Loading structural data...\n');
STRUCTURE = zeros(Nsub,1,Nvoxels);
for i = 1:Nsub
    V = spm_vol(STRUCTUREPATHS{i});
    I = spm_read_vols(V);
    STRUCTURE(i,1,:) = I(Indices);
end

% Behavioral Data
BEHAVIOR = ECFDualCORmedRT - ECFSingCORmedRT;
Thresh = [0.05 0.01 0.005];
fprintf(1,'Done prepapring data.\n');

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
%Inter([1 3],4) = 1;
Paths(1,2,1) = 1;
Paths(2,3,1) = 2;
Paths(3,4,1) = 3;

% Model 4 Parallel
Direct = zeros(ModelInfo.Nvar);
Inter = zeros(ModelInfo.Nvar);
Paths = zeros(ModelInfo.Nvar,ModelInfo.Nvar,2);
Direct(1,[2 3 4]) = 1;
Direct(2,4) = 1;
Direct(3,4) = 1;
Paths(1,2,1) = 1;
Paths(2,4,1) = 2;
Paths(1,3,2) = 1;
Paths(3,4,2) = 2;



Nboot = 0;
Nperm = 1000;
Thresh = [0.05 0.01 0.001];

ModelInfo = {};
ModelInfo.Names = Names;
ModelInfo.data = data;
ModelInfo.Direct = Direct;
ModelInfo.Inter = Inter;
ModelInfo.Paths = Paths;
ModelInfo.Nboot = Nboot;
ModelInfo.Nperm = Nperm;
ModelInfo.BaseDir = BaseDir;
ModelInfo.Tag= 'JASON';
ModelInfo.Indices = Indices;
ModelInfo.NJobSplit = 500;
ModelInfo.Thresholds = Thresh;
ModelInfo.STRAT = [];
ModelInfo.NSub = size(data{1},1);
ModelInfo.Nvar = Nvar;




%% The following is what I am preparing as the main computational piece of
% code that will be sent out to the cluster nodes.
Nvox = length(data.Indices);
Nboot = data.Nboot;
tempData = data;
Parameters = cell(Nvox,1);
Thresholds = data.Thresholds;
NThresh = length(data.Thresholds);
FieldNames = {'beta' 'B' 'Paths'};
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
    % create the structure to hold the bootstrap resample results
    for k = 1:length(FieldNames)
        Value = getfield(PointEstimateResults,FieldNames{k});
        if iscell(Value)
            BlankValue = cell(size(Value,1),Nboot);
        else
            BlankValue = zeros([size(Value) Nboot]);
        end
        BootStrap = setfield(BootStrap,FieldNames{k},BlankValue);
    end
    % run the boot strapping
    for j = 1:Nboot
        BStempData = tempData;
        BStempData.data = tempData.data(data.BootSamp(:,j),:);
        BSResults = WIPsubfnFitModel(BStempData);
        BootStrap.beta(:,:,j) = BSResults.beta;
        BootStrap.B(:,:,j) = BSResults.B;
        BootStrap.Paths{j} = BSResults.Paths;
    end
    % Calculate the statistics on the BS results.
    % I am thinking it might be a good idea to save the BCa parameters like
    % the theta value from the jack-knife step.
    % Jack-knife
    Value = getfield(PointEstimateResults,'Paths');
    JackKnifeValue = cell(size(Value,1),NSub);
    JKtempData = tempData;
    for k = 1:NSub
        Include = ones(NSub,1);
        Include(k) = 0;
        JKtempData.data = tempData.data(logical(Include),:);
        JKResults = WIPsubfnFitModel(JKtempData);
        JackKnifeValue{k} = JKResults.Paths;
    end
    NPaths = size(PointEstimateResults.Paths,1);

    [M N P] = size(PointEstimateResults.Paths{1});
    Alpha1 = zeros(NPaths,M,N,P,NThresh);
    Alpha2 = zeros(NPaths,M,N,P,NThresh);
    Zval = zeros(NPaths,M,N,P);
    pval = zeros(NPaths,M,N,P);
    
    for j = 1:NPaths
        % here I need to cycle over all values in the path. If there is a
        % two-way interaction in the path ( a single moderator) then I need to
        % cycle over a vector. If there is a three-way interaction then there
        % is a plane to step over. If there is a four-way interaction then
        % there is a 3-D space to cycle over. The dimensions are determined by
        % the size of the paths.
        
        for m = 1:M
            for n = 1:N
                for p = 1:P
                    % Extract the botstrap resamples for this path and step
                    bstat = zeros(Nboot,1);
                    for k = 1:Nboot
                        bstat(k) = BootStrap.Paths{k}{j}(m,n,p);
                    end
                    % Extract the jack-knife results for this path and step
                    theta = zeros(NSub,1);
                    for k = 1:NSub
                        theta(k) = JackKnifeValue{k}{j}(m,n,p);
                    end
                    ThetaDiff = sum(theta)/NSub - theta;
                    % find the acceleration factor
                    acc = (sum(ThetaDiff.^3))/(6*(sum(ThetaDiff.^2))^(3/2));
                    zh0 = norminv(length(find(bstat < PointEstimateResults.Paths{j}(m,n,p)))/Nboot);
                    % find the confidence intervals for the chosen thresholds
                    for t = 1:NThresh
                        zA = norminv(Thresholds(t)/2);
                        z1mA = norminv(1 - Thresholds(t)/2);
                        Alpha1(j,m,n,p,t) = normcdf(zh0 + (zh0+zA)/(1 - acc*(zh0 + zA)));
                        Alpha2(j,m,n,p,t) = normcdf(zh0 + (zh0+z1mA)/(1 - acc*(zh0 + z1mA)));
                    end
                    % Find percentile of the distribution below the null value
                    PCTlower = sum(bstat < PointEstimateResults.Paths{j}(m,n,p))./Nboot;
                    PCTupper = sum(bstat > PointEstimateResults.Paths{j}(m,n,p))./Nboot;
                    IsLower = PCTlower < PCTupper;
                    
                    
                    % Check to make sure the tails are not zero
                    PCTlower = max(PCTlower,1./Nboot);
                    PCTupper = min(PCTupper, 1-1./Nboot);
                    
                    PCTupper = max(PCTupper,1/Nboot);
                    PCTlower = min(PCTlower,1-1/Nboot);
                    
                    % Z-score
                    ZPCTlower = norminv(PCTlower) - zh0;
                    ZPCTupper = norminv(PCTupper) - zh0;
                    % adjust for acceleration
                    Zlower = (ZPCTlower.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTlower);
                    Zupper = (ZPCTupper.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTupper);
                    % which tail are we in?
                    
                    ZTemp = Zupper;
                    ZTemp(IsLower) = -Zlower(IsLower);
                    Zval(j,m,n,p) = ZTemp;
                    
                    pTemp = normcdf(ZTemp);
                    pTemp = [pTemp; 1 - pTemp];
                    pval(j,m,n,p) = 2.*min(pTemp);
                end
            end
        end
    end

    Parameters{i}.Zval = Zval;
    Parameters{i}.pval = pval;
    Parameters{i}.Alpha1 = Alpha1;
    Parameters{i}.Alpha2 = Alpha2;
    fprintf(1,'Voxel %d of %d in %0.2f seconds.\n',i,Nvox,toc);
end
%%
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

function WIPProcessBatch(BaseDir,count,Nboot)
%BaseDir = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/Model75_XAgeGroup_MpGM_MfMRI_YPerf_COVNboot0_25-Feb-2014_11-46';
Nboot = str2num(Nboot);
cd(BaseDir)
OutputFolder = fullfile(BaseDir,'BootStraps');
if ~exist(OutputFolder)
    mkdir(OutputFolder);
end
load AllData
Nvox = size(AllData.Indices,1);
load PointEstimateResults
PointEstimate = PointEstResults{1};
FieldNames = {'beta' 'B' 'Paths'};
BootStrap = {};

for i = 1:length(FieldNames)
    Value = getfield(PointEstimate,FieldNames{i});
    if iscell(Value)
        BlankValue = cell(size(Value),Nvox,Nboot);
    else
        BlankValue = zeros([size(Value) Nvox Nboot]);
    end
    BootStrap = setfield(BootStrap,FieldNames{i},BlankValue);
end
NPaths = size(PointEstimate.Paths,1);
for i = 1:Nboot
    tic
    
    BootParameters = WIPsubfnVoxelWiseProcessBatch(AllData,1);
    for j = 1:Nvox
        BootStrap.beta(:,:,j,i) = BootParameters{j}.beta;
        BootStrap.B(:,:,j,i) = BootParameters{j}.B;
        for k = 1:NPaths
            BootStrap.Paths{k,j,i} = BootParameters{j}.Paths{k};
        end
    end
    t = toc;
    fprintf(1,'Finished bootstrap %d of %d in %0.3f seconds.\n',i,Nboot,t);
end

OutFile = fullfile(OutputFolder,sprintf('BootStrapResults_%s_N%d',count,Nboot));
Str = sprintf('save %s BootStrap',OutFile);
eval(Str)
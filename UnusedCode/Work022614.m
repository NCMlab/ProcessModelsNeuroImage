SelectedDir = spm_select(1,'dir');
%SelectedDir = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/Model75_XAgeGroup_MpGM_MfMRI_YPerf_COVNboot0_25-Feb-2014_11-46';
cd(SelectedDir)
F = dir('BootStrap*.mat');
Njobs = length(F);
i = 1;
fprintf(1,'Working on chunk %04d of %d\n',i,Njobs);
load(F(i).name)
NbootSplit = size(BootStrap.beta,4);
Nvox = size(BootStrap.beta,3);
Nboot = NbootSplit*Njobs;
BSPaths = zeros(size(BootStrap.Paths{1,1,1}),Nboot,Nvox);
% extract the first chunk
for j = 1:NbootSplit
    for k = 1:Nvox
        BSPaths(:,j,k) = BootStrap.Paths{:,k,j};
    end
end

for i = 2:Njobs
    fprintf(1,'Working on chunk %04d of %d\n',i,Njobs);
    load(F(i).name)
    for j = 1:NbootSplit
        for k = 1:Nvox
            BSPaths(:,(i-1)*NbootSplit+j,k) = BootStrap.Paths{:,k,j};
        end
    end
    clear BootStrap
end

load PointEstimateResults
PEPaths = zeros([size(PointEstResults{1}.Paths{1}) Nvox]);
for k = 1:Nvox
    PEPaths(:,:,k) = PointEstResults{k}.Paths{1};
end
PEPath1 = squeeze(PEPaths(1,1,:));
PEPath2 = squeeze(PEPaths(2,1,:));
Path1 = squeeze(BSPaths(1,:,:));
Path2 = squeeze(BSPaths(2,:,:));

%%
minPath1 = min((Path1));
maxPath1 = max((Path1));

LargestValuePath1 = zeros(Nboot,1);
for i = 1:Nboot
    minV = min(Path1(i,:));
    maxV = max(Path1(i,:));
    if abs(maxV) > abs(minV)
        LargestValuePath1(i) = maxV;
    else
        LargestValuePath1(i) = minV;
    end
end
%%    
close all

maxPath1 = min((Path1));
SmaxPath1 = sort(maxPath1,2,'ascend');
maxPath2 = min((Path2));
SmaxPath2 = sort(maxPath2,2,'ascend');
a = 0.025;
c = floor(a*Nboot);
Limit1 = SmaxPath1(c)
Limit2 = SmaxPath2(c)


figure(1)
clf
hist(maxPath1,50)
h = line([Limit1 Limit1],[0 500]);
set(h,'Color','r')
figure(2)
clf
hist(PEPath1,50);
h = line([Limit1 Limit1],[0 800]);
set(h,'Color','r')

figure(3)
clf
hist(maxPath2,50)
h = line([Limit2 Limit2],[0 500]);
set(h,'Color','r')
figure(4)
clf
hist(PEPath2,50);
h = line([Limit2 Limit2],[0 800]);
set(h,'Color','r')


F1 = find(PEPath1>Limit1);
length(F1)
F2 = find(PEPath2>Limit2);
length(F2)
%%
Vo = AnalysisParameters.V;
Vo.fname = 'Path1_signMaxVoxel.nii'
I1 = zeros(Vo.dim);
I1(AnalysisParameters.Indices(F1)) = 1;
spm_write_vol(Vo,I1)

Vo = AnalysisParameters.V;
Vo.fname = 'Path2_signMaxVoxel.nii'
I1 = zeros(Vo.dim);
I1(AnalysisParameters.Indices(F2)) = 1;
spm_write_vol(Vo,I1)

Vo = AnalysisParameters.V;
Vo.fname = 'mask.nii'
I1 = zeros(Vo.dim);
I1(AnalysisParameters.Indices) = 1;
spm_write_vol(Vo,I1)


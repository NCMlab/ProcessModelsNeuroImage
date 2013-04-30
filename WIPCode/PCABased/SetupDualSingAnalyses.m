BaseDir = '/Users/jason/Documents/MyData/ECFTaskMed_PCA/Data';
SingDir = 'Sing_XageGr_MfMRI_YmedRT';
DualDir = 'Dual_XageGr_MfMRI_YmedRT';

SingData = load(fullfile(BaseDir,SingDir,'Data.mat'));
DualData = load(fullfile(BaseDir,DualDir,'Data.mat'));

NSub = length(SingData.AllData.X);

PerfDual = DualData.AllData.Y;
PerfSing = SingData.AllData.Y;

fMRIDual = squeeze(DualData.AllData.M);
fMRISing = squeeze(SingData.AllData.M);

AgeGroup = SingData.AllData.X;
Indices = SingData.AllData.Indices;
DIM = SingData.AllData.DIM;
names = SingData.AllData.names;
Vi = SingData.AnalysisParameters.V;
NVox = length(Indices);
NPCs = 12;

clear SingData;
clear DualData;

% Create matrix of all possible cominations of PCs
combo_matrix = boolean_enumeration_f(NPCs);
% how many combos are there?
NCombos = size(combo_matrix,1);
remove_row_means = 1;
% Combine the fMRI data from both Sing and Dual
fMRI_SingDual = zeros(2*NSub,NVox);
for i = 1:NSub
    fMRI_SingDual((i-1)*2+1,:) = fMRISing(i,:);
    fMRI_SingDual(i*2,:) = fMRIDual(i,:);
end

clear fMRISing
clear fMRIDual

% calculate the PCA 

[lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(fMRI_SingDual)', remove_row_means);
ssf = squeeze(fMRI_SingDual) * eigenimages_noZeroes;
ssfSubset = ssf(:,1:NPCs);


% create the subject effect design matrix
SubEffectDesign = zeros(NSub*2,NSub);
for i = 1:NSub
    SubEffectDesign((i-1)*2+1:i*2,i) = 1;
end
pInvSubEffect = pinv(SubEffectDesign);

% remove subject effects
ssfSubsetNoSubEffects = zeros(size(ssfSubset));
for i = 1:NPCs
    beta = pInvSubEffect*ssfSubset(:,i);
    fit = SubEffectDesign*beta;
    ssfSubsetNoSubEffects(:,i) = ssfSubset(:,i) - fit;
end


figure(2)
imagesc(SubEffectDesign)

figure(1)
bar(lambdas)







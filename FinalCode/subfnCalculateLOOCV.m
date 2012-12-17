function subfnCalculateLOOCV(InDataFile,NPCs)
%
NPCs = str2num(NPCs);
[PathName FileName] = fileparts(InDataFile);
OutDir = fullfile(PathName,'jobs');
% Find out what subject to leave out based on the job ID
[a b] = unix('echo $SGE_TASK_ID');
Subject = str2num(b);
fprintf(1,'Working on Subject %d\n',Subject);

load(InDataFile)
[NSub NMed NVox] = size(data.M);
NCOV = size(data.COV,2);
CurrentSubjects = [1:NSub];
CurrentSubjects(Subject) = 0;
CurrentSubjects = find(CurrentSubjects);
TrainData = [];
TrainData.X = data.X(CurrentSubjects);
TrainData.M = squeeze(data.M(CurrentSubjects,:,:));
TrainData.Y = data.Y(CurrentSubjects);
if NCOV
    TrainData.COV = data.COV(CurrentSubjects);
else
    TrainData.COV = [];
end
TrainData.ModelNum = data.ModelNum;

TestData = [];
TestData.X = data.X(Subject);
TestData.M = squeeze(data.M(Subject,:,:));
TestData.Y = data.Y(Subject);
if NCOV
    TestData.COV = data.COV(Subject,:);
end
TestData.ModelNum = data.ModelNum;

combo_matrix = boolean_enumeration_f(NPCs);
NCombos = size(combo_matrix,1);
remove_row_means = 1;
% apply PCA
% but squeeze out the multiple mediator dimension
[lambdas, eigenimages_noZeroes, w] = pca_f(TrainData.M', remove_row_means);
TrainSSF = squeeze(TrainData.M) * eigenimages_noZeroes;
TrainSSFsubset = TrainSSF(:,1:NPCs);
% cycle over all combinations of PCs

beta3 = subfnregress(TrainData.Y,[TrainData.X TrainData.COV]);
c = beta3(end - NCOV);
LOOmediationMatrix = zeros(1,NCombos);
for j = 1:NCombos
    selected_PCs = find(combo_matrix(j,:));
    NPCsSubset = length(selected_PCs);
    testdata = TrainData;
    testdata.M = TrainSSFsubset(:,selected_PCs);
    
    beta2 = subfnCallRegressPCs(testdata);
    cP = beta2(end - NCOV);
    ab = c - cP;
    %LOOerrorMatrix(1, j) = (TestData.Y...
        %- squeeze(TestData.M)'*eigenimages_noZeroes(:,selected_PCs)*b_subset(2:1+length(selected_PCs))...
        %- TestData.X*b_subset(1))^2;
    LOOmediationMatrix (1,j) = ab;
end
outFile = fullfile(OutDir,['LOOCV_' num2str(Subject)])
Str = sprintf('save %s %s' , outFile, 'LOOmediationMatrix');
eval(Str);
fprintf(1,'Saved file: %s\n',outFile);

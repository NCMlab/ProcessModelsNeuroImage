function [OutData, index] = subfnCreateOutDataStructureForModels(AllParameters,AnalysisParameters)
index = 1;
Nmed = AnalysisParameters.Nmed;
Nvoxels = length(AllParameters);
OutData = {};

% Find the first non-empty parameter set
count = 1;
while isempty(AllParameters{count})
    count = count + 1;
end

for j = 1:Nmed
    % Model 1
    Model1FieldNames = fieldnames(AllParameters{count}.Model1{j});

    for k = 1:length(Model1FieldNames)
        if isempty(strmatch(Model1FieldNames(k),'Outcome')) && isempty(strmatch(Model1FieldNames(k),'Model'))
            OutData{index}.name = ['Model1_Med' num2str(j) '_' Model1FieldNames{k} '_beta'];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['Model1{' num2str(j) '}.' Model1FieldNames{k} '.beta'];
            OutData{index}.dataType = 16;
            index = index + 1;
            OutData{index}.name = ['Model1_Med' num2str(j) '_' Model1FieldNames{k} '_t'];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['Model1{' num2str(j) '}.' Model1FieldNames{k} '.t'];
            OutData{index}.dataType = 16;
            index = index + 1;
        end
    end
end
% Model 2
Model2FieldNames = fieldnames(AllParameters{count}.Model2);

for k = 1:length(Model2FieldNames)
    if isempty(strmatch(Model2FieldNames(k),'Outcome')) && isempty(strmatch(Model2FieldNames(k),'Model'))
        OutData{index}.name = ['Model2_' Model2FieldNames{k} '_beta'];
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model2.' Model2FieldNames{k} '.beta'];
        OutData{index}.dataType = 16;
        index = index + 1;
        OutData{index}.name = ['Model2_' Model2FieldNames{k} '_t'];
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model2.' Model2FieldNames{k} '.t'];
        OutData{index}.dataType = 16;
        index = index + 1;
    end
end
% Model 3
Model3FieldNames = fieldnames(AllParameters{count}.Model3);
for k = 1:length(Model3FieldNames)
    if isempty(strmatch(Model3FieldNames(k),'Outcome')) && isempty(strmatch(Model3FieldNames(k),'Model'))
        OutData{index}.name = ['Model3_' Model3FieldNames{k} '_beta'];
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model3.' Model3FieldNames{k} '.beta'];
        OutData{index}.dataType = 16;
        index = index + 1;
        OutData{index}.name = ['Model3_' Model3FieldNames{k} '_t'];
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model3.' Model3FieldNames{k} '.t'];
        OutData{index}.dataType = 16;
        index = index + 1;
    end
end

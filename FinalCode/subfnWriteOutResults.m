function [tVoxelIndices tImageVoxelIndices] = subfnWriteOutResults(AllParameters,AnalysisParameters,OutputPath)

V = AnalysisParameters.V;
V.fname = OutputPath;
ModelNum = AnalysisParameters.ModelNum;
Tag = AnalysisParameters.Tag;

switch ModelNum
    case '4'
        OutName = [Tag '_Model' ModelNum '_'];
        count = 1;
        while isempty(AllParameters{count})
            count = count + 1;
        end
        Nvoxels = length(AllParameters);
        Nmed = AnalysisParameters.Nmed;
        
        Thresholds = AnalysisParameters.Thresholds;
        Nthr = length(Thresholds);
        VoxelIndices = zeros(Nvoxels,1);
        ImageVoxelIndices = zeros(Nvoxels,1);
        [OutData index] = subfnCreateOutDataStructureForModels(AllParameters);
        
        % Conditional Effects
        for j = 1:Nmed
           for i = 1:Nthr
                thrStr = num2str(Thresholds(i));
                OutData{index}.name = ['ABMed' num2str(j) 'sign_' num2str(Thresholds(i))];
                OutData{index}.data = zeros(Nvoxels,1);
                OutData{index}.field = ['AB' num2str(j) '{' num2str(j) '}.BCaci.alpha' thrStr(3:end)];
                OutData{index}.dataType = 2;
                index = index + 1;
            end
        end
        OutData{index}.name = 'k2';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['k2.pointEst'];
    case '7'
        
        OutName = [Tag '_Model' ModelNum '_'];
        count = 1;
        while isempty(AllParameters{count})
            count = count + 1;
        end
        Nvoxels = length(AllParameters);
        Nmed = AnalysisParameters.Nmed;
        
        Thresholds = AnalysisParameters.Thresholds;
        Nthr = length(Thresholds);
        VoxelIndices = zeros(Nvoxels,1);
        ImageVoxelIndices = zeros(Nvoxels,1);
        index = 1;
        OutData = {};
        
        for j = 1:Nmed
            % get A model variable name
            OutData{index}.name = ['AeffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.beta'];
            OutData{index}.field = ['Model1{' num2str(j) '}.X.beta'];
            index = index + 1;
            
            OutData{index}.name = ['AtMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.t'];
            OutData{index}.field = ['Model1{' num2str(j) '}.X.t'];
            index = index + 1;
            
            OutData{index}.name = ['ApMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.p'];
            OutData{index}.field = ['Model1{' num2str(j) '}.X.p'];
            index = index + 1;
            % Write out the modulator effect
            OutData{index}.name = ['WeffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.beta'];
            OutData{index}.field = ['Model1{' num2str(j) '}.W.beta'];
            index = index + 1;
            
            OutData{index}.name = ['WtMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.t'];
            OutData{index}.field = ['Model1{' num2str(j) '}.W.t'];
            index = index + 1;
            
            OutData{index}.name = ['WpMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.p'];
            OutData{index}.field = ['Model1{' num2str(j) '}.W.p'];
            index = index + 1;
            % Write out the interaction effect
            OutData{index}.name = ['XxWeffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.beta'];
            OutData{index}.field = ['Model1{' num2str(j) '}.X_x_W.beta'];
            index = index + 1;
            
            OutData{index}.name = ['XxWtMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.t'];
            OutData{index}.field = ['Model1{' num2str(j) '}.X_x_W.t'];
            index = index + 1;
            
            OutData{index}.name = ['XxWpMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['A{' num2str(j) '}.p'];
            OutData{index}.field = ['Model1{' num2str(j) '}.X_x_W.p'];
            index = index + 1;
            
            OutData{index}.name = ['BeffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['B{' num2str(j) '}.beta'];
            OutData{index}.field = ['Model2.M.beta'];
            index = index + 1;
            
            OutData{index}.name = ['BtMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['B{' num2str(j) '}.t'];
            OutData{index}.field = ['Model2.M.t'];
            index = index + 1;
            
            OutData{index}.name = ['BpMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['B{' num2str(j) '}.p'];
            OutData{index}.field = ['Model2.M.p'];
            index = index + 1;
            
            OutData{index}.name = ['CondABeffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            %OutData{index}.field = ['AB{' num2str(j) '}.pointEst'];
            OutData{index}.field =  ['CondAB1{' num2str(j) '}.pointEst'];
            index = index + 1;
            
            OutData{index}.name = ['CondABpValue' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['CondAB1{' num2str(j) '}.pValue'];
            index = index + 1;
            
            for i = 1:Nthr
                OutData{index}.name = ['CondABMed' num2str(j) 'sign_' num2str(Thresholds(i))];
                OutData{index}.data = zeros(Nvoxels,1);
                %OutData{index}.field = ['AB{' num2str(j) '}.BCaci.' num2str(Thresholds(i))];
                thrStr = num2str(Thresholds(i));
                OutData{index}.field = ['CondAB1{' num2str(j) '}.BCaci.alpha' thrStr(3:end)];
                index = index + 1;
            end
            
        end
        OutData{index}.name = 'Ceff';
        OutData{index}.data = zeros(Nvoxels,1);
        %OutData{index}.field = ['C.beta'];
        %OutData{index}.field = ['Model3.X.beta'];
        OutData{index}.field = ['Model3.X.beta'];
        index = index + 1;
        
        OutData{index}.name = 'Ct';
        OutData{index}.data = zeros(Nvoxels,1);
        %OutData{index}.field = ['C.t'];
        OutData{index}.field = ['Model3.X.t'];
        index = index + 1;
        
        OutData{index}.name = 'Cp';
        OutData{index}.data = zeros(Nvoxels,1);
        %OutData{index}.field = ['C.p'];
        OutData{index}.field = ['Model3.X.p'];
        index = index + 1;
        
        OutData{index}.name = 'CPeff';
        OutData{index}.data = zeros(Nvoxels,1);
        %OutData{index}.field = ['CP.beta'];
        OutData{index}.field = ['Model2.X.beta'];
        index = index + 1;
        
        OutData{index}.name = 'CPt';
        OutData{index}.data = zeros(Nvoxels,1);
        %OutData{index}.field = ['CP.t'];
        OutData{index}.field = ['Model2.X.t'];
        index = index + 1;
        
        OutData{index}.name = 'CPp';
        OutData{index}.data = zeros(Nvoxels,1);
        %OutData{index}.field = ['CP.p'];
        OutData{index}.field = ['Model2.X.p'];
        index = index + 1;
        
        
    case '14'
        OutName = [Tag '_Model' ModelNum '_'];
        count = 1;
        while isempty(AllParameters{count})
            count = count + 1;
        end
        Nvoxels = length(AllParameters);
        Nmed = size(AllParameters{count}.AB,2);
        Thresholds = fieldnames(AllParameters{count}.AB{1}.BCaci);
        Nthr = length(Thresholds);
        VoxelIndices = zeros(Nvoxels,1);
        ImageVoxelIndices = zeros(Nvoxels,1);
        index = 1;
        OutData = {};
        
        for j = 1:Nmed
            OutData{index}.name = ['InteffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['Int{' num2str(j) '}.beta'];
            index = index + 1;
            OutData{index}.name = ['InttMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['Int{' num2str(j) '}.t'];
            index = index + 1;
            OutData{index}.name = ['IntpMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['Int{' num2str(j) '}.p'];
            index = index + 1;
            
            OutData{index}.name = ['AeffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['A{' num2str(j) '}.beta'];
            index = index + 1;
            OutData{index}.name = ['AtMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['A{' num2str(j) '}.t'];
            index = index + 1;
            OutData{index}.name = ['ApMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['A{' num2str(j) '}.p'];
            index = index + 1;
            OutData{index}.name = ['BeffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['B{' num2str(j) '}.beta'];
            index = index + 1;
            OutData{index}.name = ['BtMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['B{' num2str(j) '}.t'];
            index = index + 1;
            OutData{index}.name = ['BpMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['B{' num2str(j) '}.p'];
            index = index + 1;
            OutData{index}.name = ['ABeffMed' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['AB{' num2str(j) '}.pointEst'];
            index = index + 1;
            OutData{index}.name = ['ABpValue' num2str(j)];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['AB{' num2str(j) '}.pValue'];
            index = index + 1;
            
            for i = 1:Nthr
                OutData{index}.name = ['ABMed' num2str(j) 'sign_' Thresholds{i}];
                OutData{index}.data = zeros(Nvoxels,1);
                OutData{index}.field = ['AB{' num2str(j) '}.BCaci.' Thresholds{i}];
                index = index + 1;
            end
        end
        OutData{index}.name = 'Ceff';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['C.beta'];
        index = index + 1;
        OutData{index}.name = 'Ct';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['C.t'];
        index = index + 1;
        OutData{index}.name = 'Cp';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['C.p'];
        index = index + 1;
        OutData{index}.name = 'CPeff';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['CP.beta'];
        index = index + 1;
        OutData{index}.name = 'CPt';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['CP.t'];
        index = index + 1;
        OutData{index}.name = 'CPp';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['CP.p'];
        index = index + 1;
        OutData{index}.name = 'k2';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['k2.pointEst'];
%% === MODEL 74 ==========================================
    case '74'
        OutName = [Tag '_Model' ModelNum '_'];
        count = 1;
        while isempty(AllParameters{count})
            count = count + 1;
        end
        Nvoxels = length(AllParameters);
        Nmed = AnalysisParameters.Nmed;
        
        Thresholds = AnalysisParameters.Thresholds;
        Nthr = length(Thresholds);
        VoxelIndices = zeros(Nvoxels,1);
        ImageVoxelIndices = zeros(Nvoxels,1);

        [OutData index] = subfnCreateOutDataStructureForModels(AllParameters,AnalysisParameters);
        % Conditional Effects
        for j = 1:Nmed
           for i = 1:Nthr
                thrStr = num2str(Thresholds(i));
                OutData{index}.name = ['CondABMed' num2str(j) 'sign_' num2str(Thresholds(i))];
                OutData{index}.data = zeros(Nvoxels,1);
                OutData{index}.field = ['CondAB' num2str(j) '{' num2str(j) '}.BCaci.alpha' thrStr(3:end)];
                OutData{index}.dataType = 2;
                index = index + 1;
            end
        end
end

% put the data into the output structure
for i = 1:Nvoxels
    if ~isempty(AllParameters{i})
        VoxelIndices(i) = 1;
        ImageVoxelIndices(i) = AnalysisParameters.Indices(i);
        for j = 1:length(OutData)
            % the confidence intervals need to be treated special
            if strfind(OutData{j}.field,'alpha')
                thresholds = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                if prod(thresholds)>0
                    OutData{j}.data(i) = 1;
                    OutData{j}.dataType = 2;
                end
            else
                % fprintf(1,'%s\n', ['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                OutData{j}.data(i) = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
            end
        end
    end
end

% write images to file
tVoxelIndices = find(VoxelIndices);
tImageVoxelIndices = ImageVoxelIndices(find(ImageVoxelIndices));
for i = 1:length(OutData)
    Vo = V;
    Vo.fname = fullfile(OutputPath,[OutName OutData{i}.name '.nii']);
    Vo.descrip = '';
    Vo.n = [1 1];
    Vo.dt = [OutData{i}.dataType 0];
    Y = zeros(Vo.dim);
    Y(tImageVoxelIndices) = OutData{i}.data(tVoxelIndices);
    spm_write_vol(Vo,Y);
end
VoxelsProcessed =  length(tVoxelIndices)


function [OutData index] = subfnCreateOutDataStructureForModels(AllParameters,AnalysisParameters)
index = 1;
Nmed = AnalysisParameters.Nmed;
Nvoxels = length(AllParameters);
OutData = {};
for j = 1:Nmed
    % Model 1
    Model1FieldNames = fieldnames(AllParameters{1}.Model1{j});
    for k = 1:length(Model1FieldNames)
        if isempty(strmatch(Model1FieldNames(k),'Outcome')) && isempty(strmatch(Model1FieldNames(k),'Model'))
            OutData{index}.name = ['Model1_Med' num2str(j) '_' Model1FieldNames{k} '_beta'];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['Model1{' num2str(j) '}.' Model1FieldNames{k} '.beta'];
            OutData{index}.dataType = 16;
            index = index + 1;
            OutData{index}.name = ['Model1_Med' num2str(j) '_' Model1FieldNames{k} '_t'];
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['Model1{' num2str(j) '}.' Model1FieldNames{k} '.t']
            OutData{index}.dataType = 16;
            index = index + 1;
        end
    end
end
% Model 2
Model2FieldNames = fieldnames(AllParameters{1}.Model2);
for k = 1:length(Model2FieldNames)
    if isempty(strmatch(Model2FieldNames(k),'Outcome')) && isempty(strmatch(Model2FieldNames(k),'Model'))
        OutData{index}.name = ['Model2_' Model2FieldNames{k} '_beta'];
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model2.' Model2FieldNames{k} '.beta'];
        OutData{index}.dataType = 16;
        index = index + 1;
        OutData{index}.name = ['Model2_' Model2FieldNames{k} '_t'];
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model2.' Model2FieldNames{k} '.t']
        OutData{index}.dataType = 16;
        index = index + 1;
    end
end
% Model 3
Model3FieldNames = fieldnames(AllParameters{1}.Model3);
for k = 1:length(Model3FieldNames)
    if isempty(strmatch(Model3FieldNames(k),'Outcome')) && isempty(strmatch(Model3FieldNames(k),'Model'))
        OutData{index}.name = ['Model3_' Model3FieldNames{k} '_beta'];
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model3.' Model3FieldNames{k} '.beta'];
        OutData{index}.dataType = 16;
        index = index + 1;
        OutData{index}.name = ['Model3_' Model3FieldNames{k} '_t'];
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model3.' Model3FieldNames{k} '.t']
        OutData{index}.dataType = 16;
        index = index + 1;
    end
end

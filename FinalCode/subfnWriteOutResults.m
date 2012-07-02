function [tVoxelIndices tImageVoxelIndices] = subfnWriteOutResults(AllParameters,AnalysisParameters)

V = AnalysisParameters.V;
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
        Nmed = size(AllParameters{count}.AB,2);
        Thresholds = fieldnames(AllParameters{count}.AB{1}.BCaci);
        Nthr = length(Thresholds);
        VoxelIndices = zeros(Nvoxels,1);
        ImageVoxelIndices = zeros(Nvoxels,1);
        index = 1;
        OutData = {};

        for j = 1:Nmed
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
 

        for i = 1:Nvoxels
            if ~isempty(AllParameters{i})
                VoxelIndices(i) = 1;
                ImageVoxelIndices(i) = AllParameters{i}.VoxelIndex;
                for j = 1:length(OutData)
                    % the confidence intervals need to be treated special
                    if strfind(OutData{j}.field,'alpha')
                        thresholds = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                        
                        if prod(thresholds)>0
                            OutData{j}.data(i) = 1;
                        end
                    else

                      % fprintf(1,'%s\n', ['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                        OutData{j}.data(i) = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                    end
                end
            end
        end
                

        tVoxelIndices = find(VoxelIndices);
        tImageVoxelIndices = ImageVoxelIndices(find(ImageVoxelIndices));

        for i = 1:length(OutData)
            Vo = V;
            Vo.fname = fullfile(Vo.fname,[OutName OutData{i}.name '.nii']);
            Vo.descrip = '';
            Vo.n = [1 1];
            Y = zeros(Vo.dim);
            Y(tImageVoxelIndices) = OutData{i}.data(tVoxelIndices);
            spm_write_vol(Vo,Y);
        end
        
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
 

        for i = 1:Nvoxels
            if ~isempty(AllParameters{i})
                VoxelIndices(i) = 1;
                ImageVoxelIndices(i) = AllParameters{i}.VoxelIndex;
                for j = 1:length(OutData)
                    % the confidence intervals need to be treated special
                    if strfind(OutData{j}.field,'alpha')
                        thresholds = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                        
                        if prod(thresholds)>0
                            OutData{j}.data(i) = 1;
                        end
                    else

                      % fprintf(1,'%s\n', ['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                        OutData{j}.data(i) = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                    end
                end
                
            end
        end
                

        tVoxelIndices = find(VoxelIndices);
        tImageVoxelIndices = ImageVoxelIndices(find(ImageVoxelIndices));

        for i = 1:length(OutData)
            Vo = V;
            Vo.n = [1 1];
            Vo.fname = fullfile(Vo.fname,[OutName OutData{i}.name '.nii']);
            Vo.descrip = '';
            Y = zeros(Vo.dim);
            Y(tImageVoxelIndices) = OutData{i}.data(tVoxelIndices);
            spm_write_vol(Vo,Y);
        end
    
        
    
        
end
            
            
VoxelsProcessed =  length(tVoxelIndices)
                   
                    
        
        
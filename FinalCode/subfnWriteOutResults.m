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
        [OutData index] = subfnCreateOutDataStructureForModels(AllParameters,AnalysisParameters);
        
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
        OutData{index}.dataType = 16;
        index = index + 1;
        % write out Model 2 R2
        OutData{index}.name = 'Model2_Rsq';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model2.Model.rsquare'];
        OutData{index}.dataType = 16;
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

        [OutData index] = subfnCreateOutDataStructureForModels(AllParameters,AnalysisParameters);
        % Conditional Effects
       
        for j = 1:Nmed
            for i = 1:Nthr
                thrStr = num2str(Thresholds(i));
                OutData{index}.name = sprintf('CondABMed%d_pV000_sign%0.4f',j,Thresholds(i));
                OutData{index}.data = zeros(Nvoxels,1);
                OutData{index}.field = ['CondAB' num2str(j) '{1}.BCaci.alpha' thrStr(3:end)];
                OutData{index}.dataType = 2;
                index = index + 1;
            end
            OutData{index}.name = sprintf('CondABMed%d_pV000_pointEst',j);
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['CondAB' num2str(j) '{1}.pointEst'];
            OutData{index}.dataType = 16;
            index = index + 1;
        end
        
        
        % write out Model 2 R2
        OutData{index}.name = 'Model2_Rsq';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model2.Model.rsquare'];
        OutData{index}.dataType = 16;
       
        
              
    
        
    case '14'
       
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
        NProbe = length(AllParameters{1}.CondAB1);
        for j = 1:Nmed
            for k = 1:NProbe
                probeValue = AllParameters{1}.CondAB1{k}.probeValue;
                for i = 1:Nthr
                    thrStr = num2str(Thresholds(i));
                    OutData{index}.name = sprintf('CondABMed%d_pV%0.2f_sign%0.4f',j,probeValue,Thresholds(i));
                    OutData{index}.data = zeros(Nvoxels,1);
                    OutData{index}.field = ['CondAB' num2str(j) '{' num2str(k) '}.BCaci.alpha' thrStr(3:end)];
                    OutData{index}.dataType = 2;
                    index = index + 1;
                    
                end
            end
        end
        
        % write out Model 2 R2
        OutData{index}.name = 'Model2_Rsq';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['Model2.Model.rsquare'];
        OutData{index}.dataType = 16;
end
%%
% put the data into the output structure
OutDataForProbesFlag = 0;
for i = 1:Nvoxels
    if ~isempty(AllParameters{i})
        VoxelIndices(i) = 1;
        ImageVoxelIndices(i) = AnalysisParameters.Indices(i);
        % this limit for the length of OutData maye change due to the
        % discovery of significant probe values in the results. Therefore, %
        % it must remain as a check of the length.
        if ~OutDataForProbesFlag
            for j = 1:length(OutData)
                % If an OutData structure refers to something with a probe
                % value then there needs to be a loop here that cycles over all
                % probe values and writes out the images. This needs to be done
                % here because we cannot assume that each voxel will have the
                % same moderating variable, i.e. there may be a voxelwise
                % moderator.
                
                % only add the new structures for the probe values if they
                % have NOT been made before.
                if ~isempty(strfind(OutData{j}.name,'000'))
                    
                    % This is a probe measure.
                    % Since this is a probe measure, the OutData cell is
                    % replaced with a structure, replicating itself for each of
                    % the probeValues.
                    % How many probe values are there for this voxel?
                    Nprobe = length(AllParameters{i}.CondAB1);
                    if Nprobe > 1
                        OutDataForProbesFlag = 1;
                    end
                    probeValues = zeros(Nprobe,1);
                    for k = 1:Nprobe
                        probeValues(k) = AllParameters{i}.CondAB1{k}.probeValue;
                    end
                    tempOutData = {};
                    index = 1;
                    tempName = OutData{j}.name;
                    tempField = OutData{j}.field;
                    for k = 1:Nprobe
                        tempOutData{index}.name = strrep(tempName,'000',sprintf('%0.2f',probeValues(k)));
                        tempOutData{index}.field = strrep(tempField,'000',num2str(k));
                        tempOutData{index}.data = zeros(Nvoxels,1);
                        tempOutData{index}.dataType = 2;
                        index = index + 1;
                    end
                    % remove the OutData cell and add the new structure to
                    % the end.
                    OutData{j} = [];
                    OutData = OutData(~cellfun('isempty',OutData));
                    index = length(OutData);
                    for k = 1:Nprobe
                        OutData{index+k} = tempOutData{k};
                    end
                end
            end
        end
        % This double checking may be slow, but it is done so that if you
        % are in the first voxel with probe values that the new additions
        % ot the OutData are still cycled over.
        for j = 1:length(OutData)
            % the confidence intervals need to be treated special
            if strfind(OutData{j}.field,'alpha')
                % This try/catch block is to see if there is more then one
                % probe value. If not then there is an error and it is
                % skipped. But this assumes that we know that there are
                % saved probe values.
                try
                    thresholds = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                    if prod(thresholds)>0
                        OutData{j}.data(i) = 1;
                        OutData{j}.dataType = 2;
                    end
                catch
                end
                
            else
                % fprintf(1,'%s\n', ['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                OutData{j}.data(i) = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
            end
        end
    end
end
%%
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
% write out the mask
Vo = V;
Vo.fname = fullfile(OutputPath,['mask.nii']);
Vo.descrip = '';
vo.n = [1 1];
Vo.dt = [2 0];
Y = zeros(Vo.dim);
Y(tImageVoxelIndices) = 1;
spm_write_vol(Vo,Y);


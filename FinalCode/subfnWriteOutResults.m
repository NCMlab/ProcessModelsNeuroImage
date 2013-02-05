function [tVoxelIndices tImageVoxelIndices] = subfnWriteOutResults(AllParameters,AnalysisParameters,OutputPath)
% The aim is to break this up into pieces so that it can be modular and
% work with voxel-based analyses or Freesurfer based analyses
% 1) create a structure with the thre model results
% 2) add structures for the conditional effects for the different models
% 3) Fill in the structure with the actual data
% 4) write out images
%   I think there shoudl be a flag in the AnalysisParameters saying what
%   sort of analysis this is: 1) voxel-wise, freesurfer, PCA
%   4a) this may be voxel-wise or
%   4b) freesurfer images
%
% All models will have some things in common

V = AnalysisParameters.V;
V.fname = OutputPath;
ModelNum = AnalysisParameters.ModelNum;
Tag = AnalysisParameters.Tag;
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
% This creates the structure for the base models that are common to all
% models
[OutData index] = subfnCreateOutDataStructureForModels(AllParameters,AnalysisParameters);

% write out Model 2 R-squared
OutData{index}.name = 'Model2_Rsq';
OutData{index}.data = zeros(Nvoxels,1);
OutData{index}.field = ['Model2.Model.rsquare'];
OutData{index}.dataType = 16;

switch ModelNum
    case '4'
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
        
    case '7'
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
        
    case '58'
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
        
    case '74'
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
end
%%
OutData = subfnPutDataIntoOutputStructure(OutData,AllParameters,AnalysisParameters)
%%

subfnWriteResultsToImages(VoxelIndices,OutData,OutName,ImageVoxelIndices,V)


function subfnWriteResultsToImages(VoxelIndices,OutData,OutName,ImageVoxelIndices,V)
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
VoxelsProcessed =  length(tVoxelIndices);
% write out the mask
Vo = V;
Vo.fname = fullfile(OutputPath,['mask.nii']);
Vo.descrip = '';
vo.n = [1 1];
Vo.dt = [2 0];
Y = zeros(Vo.dim);
Y(tImageVoxelIndices) = 1;
spm_write_vol(Vo,Y);


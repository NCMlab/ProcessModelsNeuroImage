function Results = CycleOverVoxelsProcessBootstrap(ModelInfo)
% The number of voxels variable refers to the number of voxels being
% analyzed in this chunk. If the data is being split into chunks for
% analysis on a cluster then this number will differe from teh total number
% of voxels in the analysis.

Nvoxels = length(ModelInfo.Indices);
% Prepare the output structure
Results = cell(Nvoxels,1);
for i = 1:Nvoxels
    fprintf(1,'%d of %d voxels\n',i,Nvoxels);
    % Extract the data for this voxel
    OneVoxelModel = ExtractDataFromVoxel(ModelInfo,i);
    % Perform the analysis for this voxel
    Results{i} = OneVoxelProcessBootstrap(OneVoxelModel);
end



function OneVoxelModel = ExtractDataFromVoxel(ModelInfo,Index)
% Extract the data for a single voxel
OneVoxelModel = ModelInfo;
OneVoxelModel.data = zeros(ModelInfo.Nsub,ModelInfo.Nvar);
for i = 1:ModelInfo.Nvar
    if size(ModelInfo.data{i},2) > 1
        % This is voxel-wise data
        OneVoxelModel.data(:,i) = ModelInfo.data{i}(:,Index);
    else
        OneVoxelModel.data(:,i) = ModelInfo.data{i};
    end
end




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

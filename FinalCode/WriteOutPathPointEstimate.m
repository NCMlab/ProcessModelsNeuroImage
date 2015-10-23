function WriteOutPathPointEstimate(Parameters,ModelInfo)

NPaths = size(Parameters{1}.Paths,1);
temp = zeros(ModelInfo.Nvoxels,NPaths);
tempT = zeros(ModelInfo.Nvoxels,NPaths);
for j = 1:ModelInfo.Nvoxels
    for k = 1:NPaths
        temp(j,k) = Parameters{j}.Paths{k};
        tempT(j,k) = Parameters{j}.PathsTnorm{k};
    end
end
for i = 1:NPaths
    FileName = fullfile(ModelInfo.ResultsPath,sprintf('Path%d_PointEstimate.nii',i));
    I = zeros(ModelInfo.DataHeader.dim);
    I(ModelInfo.Indices) = squeeze(temp(:,i));
    % Create the header for this image
    Vo = ModelInfo.DataHeader;
    Vo.fname = FileName;
    spm_write_vol(Vo,I);
end

for i = 1:NPaths
    FileName = fullfile(ModelInfo.ResultsPath,sprintf('Path%d_Tnorm.nii',i));
    I = zeros(ModelInfo.DataHeader.dim);
    I(ModelInfo.Indices) = squeeze(tempT(:,i));
    % Create the header for this image
    Vo = ModelInfo.DataHeader;
    Vo.fname = FileName;
    spm_write_vol(Vo,I);
end
function WriteOutPathPointEstimate(Parameters,ModelInfo)
% Is this path modulated?

[NPaths] = size(Parameters{1}.Paths{1},1);
% How many steps are there in the path?
NSteps = size(Parameters{1}.Paths{1},2);
temp = zeros(ModelInfo.Nvoxels,NPaths,NSteps);
tempT = zeros(ModelInfo.Nvoxels,NPaths,NSteps);

for j = 1:ModelInfo.Nvoxels
    for k = 1:NPaths
        temp(j,k,:) = Parameters{j}.Paths{1}(k);
        tempT(j,k,:) = Parameters{j}.PathsTnorm{1}(k);
    end
end
for i = 1:NPaths
    for k = 1:NSteps
        FileName = fullfile(ModelInfo.ResultsPath,sprintf('Path%d_Step%03d_PE.nii',i,k));
        I = zeros(ModelInfo.DataHeader.dim);
        I(ModelInfo.Indices) = squeeze(temp(:,i,k));
        % Create the header for this image
        Vo = ModelInfo.DataHeader;
        Vo.fname = FileName;
        spm_write_vol(Vo,I);
    end
end

for i = 1:NPaths
    for k = 1:NSteps
        FileName = fullfile(ModelInfo.ResultsPath,sprintf('Path%d_Step%03d_Tnorm.nii',i,k));
        I = zeros(ModelInfo.DataHeader.dim);
        I(ModelInfo.Indices) = squeeze(tempT(:,i,k));
        % Create the header for this image
        Vo = ModelInfo.DataHeader;
        Vo.fname = FileName;
        spm_write_vol(Vo,I);
    end
end
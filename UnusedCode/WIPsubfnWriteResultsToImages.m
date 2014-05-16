function WIPsubfnWriteResultsToImages(OutData,ImageVoxelIndices,V,BaseName)
% write images to file

for i = 1:length(OutData)
    Vo = V;
    Vo.fname = fullfile(BaseName,[OutData{i}.name '.nii']);
    Vo.descrip = '';
    Vo.n = [1 1];
    Vo.dt = [OutData{i}.dataType 0];
    Y = zeros(Vo.dim);
    Y(ImageVoxelIndices) = OutData{i}.data;
    spm_write_vol(Vo,Y);
end

% write out the mask
Vo = V;
Vo.fname = fullfile(BaseName,['mask.nii']);
Vo.descrip = '';
Vo.n = [1 1];
Vo.dt = [2 0];
Y = zeros(Vo.dim);
Y(ImageVoxelIndices) = 1;
spm_write_vol(Vo,Y);


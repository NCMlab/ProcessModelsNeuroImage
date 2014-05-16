function subfnWriteResultsToImages(VoxelIndices,OutData,OutName,ImageVoxelIndices,V)
% write images to file
tVoxelIndices = find(ImageVoxelIndices);
tImageVoxelIndices = ImageVoxelIndices(find(ImageVoxelIndices));
for i = 1:length(OutData)
    Vo = V;
    Vo.fname = fullfile(V.fname,[OutName OutData{i}.name '.nii']);
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
Vo.fname = fullfile(V.fname,['mask.nii']);
Vo.descrip = '';
vo.n = [1 1];
Vo.dt = [2 0];
Y = zeros(Vo.dim);
Y(tImageVoxelIndices) = 1;
spm_write_vol(Vo,Y);


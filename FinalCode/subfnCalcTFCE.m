function [TFCEvalues] = subfnCalcTFCE(ImageData, ModelInfo)

% Create a temp mask image
Im = zeros(ModelInfo.DataHeader.dim);
Im(ModelInfo.Indices) = 1;
Vm = ModelInfo.DataHeader;
Vm.fname = fullfile(ModelInfo.BaseDir,'tempMask.nii');
spm_write_vol(Vm,Im);
% Create the temp stat image
I = zeros(ModelInfo.DataHeader.dim);
I(ModelInfo.Indices) = ImageData;
Vi = ModelInfo.DataHeader;
Vi.fname = fullfile(ModelInfo.BaseDir,'TestTemp.nii');
spm_write_vol(Vi,I);

% calculate TFCE on stat image
setenv('FSLOUTPUTTYPE','NIFTI');
path1 = getenv('PATH');
if isempty(strfind(path1,'fsl'))
    path1 = [path1 ':/usr/local/fsl/bin'];
    setenv('PATH', path1);
end

FSLMathsPath = 'fslmaths  ';
FSLStatsPath = 'fslstats  ';
 %[a FSLMathsPath] = unix('! which fslmaths');
% FSLMathsPath = '/usr/local/fsl/bin/fslmaths ';
 %[a FSLStatsPath] = unix('! which fslstats');
 
 
[PathName, FileName] = fileparts(Vi.fname);
outFile = fullfile(PathName,'TFCEtemp.nii');

% Positive voxels
Str = sprintf('! %s %s -mas %s -nan -tfce 2 0.5 6 %s',FSLMathsPath(1:end-1), Vi.fname,Vm.fname,outFile);
[a b] = unix(Str);
% load TFCE image
if a == 1 % Check to see if there are any errors
    Vtfce = spm_vol(outFile);
    Itfce = spm_read_vols(Vtfce);
    posTFCEvalues = Itfce(ModelInfo.Indices);
else
    posTFCEvalues = zeros(ModelInfo.Nvoxels,1);
end

% Negative voxels
% Now do it for the negative direction
Str = sprintf('! %s %s -mul -1 -thr 0 %s',FSLMathsPath(1:end-1), Vi.fname, outFile);
[a b] = unix(Str);
Str = sprintf('! %s %s -mas %s -nan -tfce 2 0.5 6 %s',FSLMathsPath(1:end-1), outFile,Vm.fname,outFile);
[a b] = unix(Str);
% load TFCE image
if a == 1 % Check to see if there are any errors
    Vtfce = spm_vol(outFile);
    Itfce = spm_read_vols(Vtfce);
    negTFCEvalues = -1.*Itfce(ModelInfo.Indices);
else
    negTFCEvalues = zeros(ModelInfo.Nvoxels,1);
end
TFCEvalues = posTFCEvalues + negTFCEvalues;


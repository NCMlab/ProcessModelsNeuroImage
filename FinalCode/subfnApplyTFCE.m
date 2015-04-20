function [maxTFCE, minTFCE] = subfnApplyTFCE(ImagePath,MaskPath)
% Take a single input image and find the max and min TFCE values using the
% FSL commands

% If only files names are given assume the file is in the current
% directory.



% Find fslmaths
setenv('FSLOUTPUTTYPE','NIFTI');
 [a FSLMathsPath] = unix('! which fslmaths');
% FSLMathsPath = '/usr/local/fsl/bin/fslmaths ';
 [a FSLStatsPath] = unix('! which fslstats');
% FSLStatsPath = '/usr/local/fsl/bin/fslstats ';

Vi = spm_vol(ImagePath);
Vm = spm_vol(MaskPath);
[PathName, FileName] = fileparts(Vi.fname);
outFile = fullfile(PathName,'TFCEtemp.nii');

Str = sprintf('! %s %s -mas %s -tfce 2 0.5 6 %s',FSLMathsPath(1:end-1), Vi.fname,Vm.fname,outFile);
[a b] = unix(Str);
if a == 1 % Check to see if there are any errors
    % Errors arise if all values in an image are negative
    % use fslstats to find the maximum in the file
    Str = sprintf('! %s %s -R',FSLStatsPath(1:end-1),outFile);
    [a b] = unix(Str);
    findUnder = find(b == ' ');
    % save this value
    posMaxTFCE = str2double(b(findUnder(1)+1:findUnder(2)-1));
    posMinTFCE = str2double(b(1:findUnder(1)-1));
else
    posMaxTFCE = 0;
    posMinTFCE = 0;
end
% Now do it for the negative direction
Str = sprintf('! %s %s -mul -1 -thr 0 %s',FSLMathsPath(1:end-1), Vi.fname, outFile);
[a b] = unix(Str);

% Use the fslmaths TFCE command on the saved file
Str = sprintf('! %s %s -mas %s -tfce 2 0.5 6 %s',FSLMathsPath(1:end-1), outFile,Vm.fname,outFile);
[a b] = unix(Str);
if a == 1    
    % use fslstats to find the maximum in the file
    Str = sprintf('! %s %s -R',FSLStatsPath(1:end-1),outFile);
    [a b] = unix(Str);
    findUnder = find(b == ' ');
    % save this value
    negMaxTFCE = str2double(b(findUnder(1)+1:findUnder(2)-1));
    negMinTFCE = str2double(b(1:findUnder(1)-1));
else
    negMaxTFCE = 0;
    negMinTFCE = 0;
end 
% Combine the two directions
if posMaxTFCE > negMinTFCE
    maxTFCE = posMaxTFCE;
else
    maxTFCE = negMinTFCE;
end

if posMinTFCE > negMaxTFCE
    minTFCE = posMinTFCE;
else
    minTFCE = -1.*negMaxTFCE;
end


function [posOutStruct negOutStruct] = subfnCreateTableOfClusterResults

Paal = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/masks/raal.nii';
Pba = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/masks/rbrodmann.nii';
% load the atlas maps
Vba = spm_vol(Pba);
Vaal = spm_vol(Paal);
Iaal = spm_read_vols(Vaal);
Iba = spm_read_vols(Vba);


% Select the thresholded image
InputImage = spm_select(1,'image');
cd(fileparts(InputImage))

choice = questdlg('Select Mask Image?', ...
	'Mask?', ...
	'yes','no','no');
switch choice
    case 'yes'
        MaskImage = spm_select(1,'image','Select mask image');
    case 'no'
        MaskImage = '';
end

choice = questdlg('Is this image thresholded?', ...
	'Thresholded?', ...
	'yes','no','yes');
switch choice
    case 'yes'
        HeightThreshold = 0.5;
        ExtentThreshold = 0.5;
    case 'no'
        prompt = {'Enter height threshold:','Enter extent threshold:'};
        dlg_title = 'Input for thresholds';
        num_lines = 1;
        def = {'1.96','100'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        HeightThreshold = str2num(answer{1});
        ExtentThreshold = str2num(answer{2});
        InputImage = subfnApplyThresholdsToImages(InputImage,HeightThreshold,ExtentThreshold,MaskImage);
end



% Create the table of results using SPM
% POSITIVE DIRECTION
[SPM xSPM] = DisplayCov(InputImage,HeightThreshold,ExtentThreshold);
NumLocalMaxima = 3;
DistancebetweenMaxima = 8;
% Check to see if all the voxelwise values are one. If so then this is an
% image of thresholded significant effects and we need to replace these
% values with effects sizes.
if unique(xSPM.Z) == 1
    EffectImage = spm_select(1,'image');
    Veffect = spm_vol(EffectImage);
    Ieffect = spm_read_vols(Veffect);
    effectZ = zeros(size(xSPM.Z));
    for i = 1:length(xSPM.Z)
        effectZ(i) = Ieffect(xSPM.XYZ(1,i),xSPM.XYZ(2,i),xSPM.XYZ(3,i));
    end
    xSPM.Z = effectZ;
end
posOutStruct = findAALandBAfromTabDat(xSPM,NumLocalMaxima,DistancebetweenMaxima,Iaal,Iba);
posOutStruct.InputImage = InputImage;
posOutStruct.VoxelsMM = xSPM.XYZmm;
posOutStruct.VoxelsVOX = xSPM.XYZ;


% NEGATIVE DIRECTION
% load the image
V = spm_vol(InputImage);
I = spm_read_vols(V);
% find the activations in the negative direction
F = find(I < -HeightThreshold);
[x y z] = ind2sub(V.dim,F);
negXYZ = [x y z];
negXYZmm = (SPM.xVol.M*[negXYZ ones(length(negXYZ),1)]')';
negZ = I(F)';
negxSPM = xSPM;
negxSPM.Z = -1.*negZ;
negxSPM.XYZ = negXYZ';
negxSPM.XYZmm = negXYZmm(:,1:3)';
negxSPM.title = 'negative direction';


negOutStruct = findAALandBAfromTabDat(negxSPM,NumLocalMaxima,DistancebetweenMaxima,Iaal,Iba);
negOutStruct.InputImage = InputImage;
negOutStruct.VoxelsMM = negxSPM.XYZmm;
negOutStruct.VoxelsVOX = negxSPM.XYZ;


fprintf(1,'===== %s ======\n','POSITIVE DIRECTION');
WriteTableOutResultsToScreen(posOutStruct,InputImage,'POS')
fprintf(1,'===== %s ======\n','NEGATIVE DIRECTION');
WriteTableOutResultsToScreen(negOutStruct,InputImage,'NEG')


%CreateMaskOfResults(posOutStruct,xSPM,negxSPM)

% Find the negative clusters

% find the locations of the maxima for each cluster

% fine the locations of the local maxima in each cluster

% find teh BA/AAL locations for these maxima
function CreateMaskOfResults(posOutStruct,xSPM,negxSPM)
[PathName FileName Ext] = fileparts(posOutStruct.InputImage);
I = zeros(xSPM.Vspm.dim);
for i = 1:length(xSPM.XYZ)
    I(xSPM.XYZ(1,i),xSPM.XYZ(2,i),xSPM.XYZ(3,i)) = 1;
end
for i = 1:length(negxSPM.XYZ)
    I(negxSPM.XYZ(1,i),negxSPM.XYZ(2,i),negxSPM.XYZ(3,i)) = 1;
end
OutFile = fullfile(PathName, [FileName '_mask' Ext]);
Vo = xSPM.Vspm;
Vo.fname = OutFile;
Vo.dt = [2 0];
spm_write_vol(Vo,I);

function WriteTableOutResultsToScreen(OutStruct,InputImage,direction)
NCl = length(OutStruct.t);
fid = 1;
fprintf(fid,'===== %s ======\n',InputImage);
fprintf(fid,'%-20s\t%5s\t%5s\t%5s\t%5s\t%5s\t%10s\t%10s\n','Region','Lat','BA','Xmm','Ymm','Zmm','Z','ClSize');
for i = 1:NCl
    fprintf(fid,'%-20s\t%5s\t%5d\t',OutStruct.aal{i},OutStruct.hemi{i}, OutStruct.ba(i));
    fprintf(fid,'%5d\t%5d\t%5d\t',OutStruct.loc(i,:));
    
     t = OutStruct.t(i);
     if strmatch(direction,'NEG')   
        fprintf(fid,'%10.3f\t', -t);
    else
        fprintf(fid,'%10.3f\t', OutStruct.t(i));
    end
        fprintf(fid,'%10s\n',OutStruct.k{i});
end
fprintf(fid,'=============================================================\n');

function [AALList BAList HemiList] = subfnLocalFindAALandBA(XYZ, Iaal, Iba)
[aalCol1 aalCol2 aalCol3] = textread('/share/studies/CogRes/GroupAnalyses/ModMedCogRes/masks/aal.nii.txt','%d%s%d');
AALlabels = FormatALLNames(aalCol2);
NVoxels = size(XYZ,1);
HemiList = {};
AALList = {};
BAList = zeros(NVoxels,1);
for i = 1:NVoxels
    CurrentValue = Iaal(XYZ(i,1), XYZ(i,2), XYZ(i,3));
    BAList(i,1) = Iba(XYZ(i,1), XYZ(i,2), XYZ(i,3));
    if CurrentValue
        %AALList{i} = aalCol2{CurrentValue};
        AALList{i} = AALlabels{CurrentValue}.out;
        HemiList{i} = AALlabels{CurrentValue}.hemi;
    else
        AALList{i} = '**empty**';
        HemiList{i} = '-';
    end
end



function OutStruct = findAALandBAfromTabDat(xSPM,NumLocalMaxima,DistancebetweenMaxima,Iaal,Iba)
TabDat = spm_list('table',xSPM,[NumLocalMaxima,DistancebetweenMaxima,'']);
% Read the table 
NClust = size(TabDat.dat,1);
XYZmm = zeros(NClust,3);
TList = zeros(NClust,1);
ZList = zeros(NClust,1);
pList = zeros(NClust,1);
ClList = {};%zeros(NClust,1);
%
Hdr = char(TabDat.hdr(2,:));
for i = 1:length(Hdr)
    element = deblank(Hdr(i,:));
    switch element
        case {'equivk'}
            KCol = i;
        case {'T'}
            TCol = i;
        case {'equivZ'}
            ZCol = i;
    end
end
for i = 1:NClust
    XYZmm(i,:) = TabDat.dat{i,end}(:)';
    TList(i,1) = TabDat.dat{i,TCol};
    ZList(i,1) = TabDat.dat{i,ZCol};
    pList(i,1) = TabDat.dat{i,end - 1};
    tempCl = TabDat.dat{i,KCol};
    if ~isempty(tempCl)
        ClList{i} = num2str(tempCl);
    else
        ClList{i} = '--';
    end
end

XYZ = (inv(xSPM.Vspm.mat)*([XYZmm ones(size(XYZmm,1),1)])')';
% This program finds the AAL and Brodmann Area labels for all XYZ mm
% coordinates it is passed.
[AALList BAList HemiList] = subfnLocalFindAALandBA(XYZ(:,1:3),Iaal,Iba);

OutStruct = {};
OutStruct.loc = XYZmm;
OutStruct.t = TList;
OutStruct.Z = ZList;
OutStruct.p = pList;
OutStruct.k = ClList';
OutStruct.aal = AALList';
OutStruct.ba = BAList;
OutStruct.hemi =  HemiList;


function AALlabels = FormatALLNames(aalCol2)
N = length(aalCol2);
AALlabels = {};
for i = 1:N
    AALlabels{i}.in = aalCol2{i};
    temp = aalCol2{i};
    fUnder = findstr(temp,'_');
    lastPiece = temp(fUnder(end)+1:end);
    if ~isempty(strmatch(lastPiece,'L')) || ~isempty(strmatch(lastPiece,'R'))
        AALlabels{i}.hemi = lastPiece;
    else
        AALlabels{i}.hemi = '';
    end
    prefix = '';
    if length(fUnder) > 1
        for j = 1:length(fUnder) - 1
            prefix = [prefix temp(fUnder(j)+1:fUnder(j+1)-1) '. '];
        end
    end
    AALlabels{i}.prefix = prefix;
    AALlabels{i}.name = temp(1:fUnder(1)-1);
   AALlabels{i}.out = [prefix AALlabels{i}.name];
end
    


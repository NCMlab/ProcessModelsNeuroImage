function [XYZmm ClList] = FindAALandBAForSPMResults(xSPM, hReg,Title, Num,Dis, OtherImage)
% Usage:
% After results have been estimated and the glass brains are displayed in
% the Graphics window. 
% Execute: 
% FindAALandBAForSPMResults(xSPM, hReg,Num,Dis,OtherImage)
% 
% NOTE: This has been written for spm5.
% 
% Extract all voxels from the SPM table of results and find the AAL and
% Brodmann Labels for them. Then write out most of the SPM results table
% and the AAL and BA labels to a tab-delimited text file.
%
% Note you will need to change the locations for MRICro which are hardcoded
% below.
%
% filename: FindAALandBAForSPMResults.m
% written by: Jason Steffener
% Date: 4/17/08

% Once SPM results are run this program will work
% The table of results is extracted

if nargin < 2
    
    fprintf('usage: FindAALandBAForSPMResults(xSPM, hReg, [Title], [Num],[Dis])\n')
    fprintf('Num and Dis are OPTIONAL arguements\n');
   fprintf('Num: number of local maxima\nDis: Distance between local maxima\n')
   
   return
elseif nargin == 2
    Title = '';
    Num = 3;
    Dis = 8;
    OtherImage = 0;
elseif nargin == 3
    
    Num = 3;
    Dis = 8;
    OtherImage = 0;
elseif nargin == 4
    Dis = 8;
    OtherImage = 0;
elseif nargin == 5
    OtherImage = 0;
end

TabDat = spm_list('List', xSPM, hReg, Num, Dis);
    %Num    - number of maxima per cluster
    %Dis    - distance among clusters (mm)
%TabDat = spm_list('ListCluster', xSPM, hReg);
% The tabular data is broken into different arrays
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

% This program finds the AAL and Brodmann Area labels for all XYZ mm
% coordinates it is passed.
[AALList BAList] = subfnFindAALandBA(XYZmm);
% Create an informative output file name
OutPutFileName = [Title xSPM.title '_Teq' num2str(xSPM.u)];
OutPutFileName(findstr(OutPutFileName, ' ')) = '';
OutPutFileName(findstr(OutPutFileName, '.')) = 'p';
OutPutFileName(findstr(OutPutFileName, '(')) = '';
OutPutFileName(findstr(OutPutFileName, ')')) = '';
OutPutFileName(findstr(OutPutFileName, '<')) = 'l';
OutPutFileName(findstr(OutPutFileName, '>')) = 'g';
OutPutFileName = [OutPutFileName '_AALandBA.csv'];
OutPutFile = fullfile(xSPM.swd,[OutPutFileName]);
fid = fopen(OutPutFile, 'w');
% Write all data to file
fprintf(fid,'%s,Hemi.,B.A.,x,y,z,Z,k\n','Region');
for i = 1:NClust
    fprintf(fid,'%s,%s,%s,',AALList{i}(1:end-2), AALList{i}(end),BAList{i});
    fprintf(fid,'%d,%d,%d,',XYZmm(i,:));
%    fprintf(fid,'%0.2f\t%0.2f\t%0.2f\t%d\n',TList(i,1),ZList(i,1),pList(i,1),ClList(i,1));
    fprintf(fid,'%0.2f,%s\n',TList(i,1),ClList{i});
end
 fclose(fid);
fprintf('Data saved to:\n\t%s\n',OutPutFile);
if OtherImage
    FindAALandBAForCurrentSPMResultsinOTHERImage(XYZmm,xSPM);
end

fprintf(1,'%20s\tLat\tBA\tx\ty\tz\tT\tk\n','Region');
for i = 1:NClust
    fprintf(1,'%20s\t%s\t%s\t',AALList{i}(1:end-2), AALList{i}(end),BAList{i});
    fprintf(1,'%d\t%d\t%d\t',XYZmm(i,:));
%    fprintf(fid,'%0.2f\t%0.2f\t%0.2f\t%d\n',TList(i,1),ZList(i,1),pList(i,1),ClList(i,1));
    fprintf(1,'%0.3f\t%s\n',TList(i,1),ClList{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AALList BAList] = subfnFindAALandBA(XYZmm)
% This takes an array of XYZ coordinates in millimeters and finds the AAL
% and BA labels for them. Using mm coordinates they are converted to the
% voxel space of teh AAL and BA maps regardless of differing image sizes.
% 
% filename: FindAALandBAForSPMResults.m
% written by: Jason Steffener
% Date: 4/17/08
%
%%%%% Change these before using!! %%%%%%%%%%%%%%%%%%%%%%%%
% Paal = '/mnt/data6/Brady/mricro/aal.img';
% [aalCol1 aalCol2 aalCol3] = textread('/mnt/data6/Brady/mricro/aal.txt','%d%s%d');
% Pba = '/mnt/data6/Brady/mricro/brodmann.img';
CodePath = fileparts(mfilename('fullpath'));

MaskPath = fullfile(fileparts(CodePath),'masks');
% Paal = '/share/data/data9/DARPA2_LS_SPM5/raal.nii';
% Pba = '/share/data/data9/DARPA2_LS_SPM5/rbrodmann.nii';
%Pcort = '/share/data/data5/locally_written_m_files/templates/rHarvardOxford-cort-prob-2mm.nii'
% if ismac
%     BaseDir = '/Users/jason/Dropbox/SteffenerColumbia/Scripts/ProcessModelsNeuroImage/masks';

    Paal = fullfile(MaskPath,'raal.nii');
    Pba = fullfile(MaskPath, 'rbrodmann.nii');
    [aalCol1 aalCol2 aalCol3] = textread(fullfile(MaskPath,'aal.nii.txt'),'%d%s%d');

% else
%     Paal = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/masks/raal.nii';
%     Pba = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/masks/rbrodmann.nii';
% 
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vba = spm_vol(Pba);
Vaal = spm_vol(Paal);
Iaal = spm_read_vols(Vaal);
Iba = spm_read_vols(Vba);
%Vcort = spm_vol(Pcort);
%Icort = spm_read_vols(Vcort);
% Account for differing image sizes and flip the images
VMat = Vaal.mat;
VMat(1,1) = -VMat(1,1);
VMat(1,4) = -VMat(1,4);
NVoxels = size(XYZmm,1);
AALList = {};
BAList = {};
for i = 1:NVoxels
    CurrentMMLoc = XYZmm(i,:)';
    CurrentAALLoc = inv(VMat)*[CurrentMMLoc; 1];
    CurrentValue = Iaal(CurrentAALLoc(1), CurrentAALLoc(2), CurrentAALLoc(3));
    BAList{i} = num2str(Iba(CurrentAALLoc(1), CurrentAALLoc(2), CurrentAALLoc(3)));
    if CurrentValue
        AALList{i} = aalCol2{CurrentValue};
    else
        if XYZmm(i,1) < 0 
            AALList{i} = '---R';
        elseif XYZmm(i,1) > 0 
            AALList{i} = '---L';
        else
            AALList{i} = '----';
        end
    end
    if ~isempty(strmatch(BAList{i},'0'))
        BAList{i} = '--';
    elseif ~isempty(strmatch(BAList{i},'48'))
        BAList{i} = '--';
    end
end



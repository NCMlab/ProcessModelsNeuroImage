function FindAALandBAFromStatImage(HeightThr, ClusterThr)
% Find the significant voxels and clusters and find the AAL and BA
% locations for them.
% filename: FindAALandBAFromStatImage.m
% date: 6/17/08
% written by: Jason Steffener
% NOTE: make sure to specify where your atlas images are
% MAKE SURE THEY ARE RESLICED TO THE CORRECT DIMENSIONS

% where are the AAL and BA maps that have been converted(resliced) into
% the data's space
Flip = inputdlg('Flip Image? (y/n)');
Paal = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/masks/raal.nii';
Pba = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/masks/rbrodmann.nii';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HERE IS A USE OF AN SPM FUNCTION
P = spm_select(1,'image','Select statistical image');
% Choose a height and extent threshold
% HeightThr = 2.3;
% ClusterThr = 40;
% HERE IS A USE OF AN SPM FUNCTION
% load up stat input image
Vin = spm_vol(P);

Img = spm_read_vols(Vin);
% find voxels above threshold
F = find(Img > HeightThr);
% convert indices to voxel locations
DIM = Vin.dim;
[x y z] = ind2sub(DIM,F);
% HERE IS A USE OF AN SPM FUNCTION
% find clusters 
[A] = spm_clusters([x y z]');
% how many clusters
NumClust = max(A);
% how big are the clusters
ClusterSize = zeros(NumClust,1);
for i = 1:NumClust
    ClusterSize(i,1) = length(find(A == i));
end
% which clusters are big enough
SignClusters = find(ClusterSize >= ClusterThr);

SignVoxelsAndClusters = zeros(size(A));
ClusterMaxStat = zeros(length(SignClusters),1);
ClusterSize = zeros(length(SignClusters),1);
ClusterMaxLoc = [];
% what is the maximal value in the clusters and how big is each cluster
for i = 1:length(SignClusters)
    F = find(A == SignClusters(i));
    tempValues = zeros(size(F));
    ClusterSize(i) = length(tempValues);
    for j = 1:length(F)
        tempValues(j) = Img(x(F(j)), y(F(j)), z(F(j)));
    end
    ClusterMaxStat(i) = max(tempValues);
    tempMaxLoc = F(find(tempValues == ClusterMaxStat(i)));
    ClusterMaxLoc(i,:) = [x(tempMaxLoc) y(tempMaxLoc) z(tempMaxLoc)];
    SignVoxelsAndClusters(F) = 1;
end
% % if strmatch(upper(Flip),'Y')
% %     warning('Flipping results');
% %     for i = 1:size(ClusterMaxLoc,1)
% %         if ClusterMaxLoc(i,1)/floor(Vin.dim(1)/2) > 1
% %             ClusterMaxLoc(i,1)  = ClusterMaxLoc(i,1) - floor(Vin.dim(1)/2);
% %         elseif ClusterMaxLoc(i,1)/floor(Vin.dim(1)/2) < 1
% %             ClusterMaxLoc(i,1)  = ClusterMaxLoc(i,1) + floor(Vin.dim(1)/2);
% % 
% %         end
% %     end
% % end
ClusterSize
% find the locatio of the cluster maximum
ClusterMaxLocmm = Vin.mat*[ClusterMaxLoc ones(size(ClusterMaxLoc,1),1)]';
% ClusterMaxLocmm = ClusterMaxLocmm(1:3,:)';
% find the location of all significant voxels within significant clusters
xSign = x(find(SignVoxelsAndClusters));
ySign = y(find(SignVoxelsAndClusters));
zSign = z(find(SignVoxelsAndClusters));

XYZ = [xSign ySign zSign];
XYZmm = Vin.mat*[XYZ ones(length(XYZ),1)]';
XYZmm = XYZmm(1:3,:)';

% load the atlas maps
Vba = spm_vol(Pba);
Vaal = spm_vol(Paal);
Iaal = spm_read_vols(Vaal);
Iba = spm_read_vols(Vba);
% find the AAL and BA locations
%[AllAALList AllBAList] = subfnLocalFindAALandBA(XYZ, Iaal, Iba);
%[MaxAALList MaxBAList] = subfnLocalFindAALandBA(XYZ, Iaal, Iba);

[MaxAALList MaxBAList] = subfnLocalFindAALandBA(ClusterMaxLoc, Iaal, Iba);

% write the data out to a file
% Create an informative output file name
[PathName FileName] = fileparts(P);
ClusterOutPutFileName = fullfile(PathName,['ClusterMax' FileName '.txt']);
fid = fopen(ClusterOutPutFileName, 'w');

% Write all data to file
if strmatch(upper(Flip),'Y')
%    warning('Flipping results');
%    Vin.mat(1,:) = -1.*Vin.mat(1,:);
fprintf(fid,'Make sure to flip Left and Right\n');
end
fprintf(fid,'%20s\t%5s\t%5s\t%5s\t%5s\t%5s\t%10s\t%10s\n','Region','Lat','BA','Xmm','Ymm','Zmm','ClMax','ClSize');
for i = 1:size(ClusterMaxLocmm,2)
    fprintf(fid,'%20s\t%5s\t%5d\t',MaxAALList{i}(1:end-2), MaxAALList{i}(end),MaxBAList(i));
    fprintf(fid,'%5d\t%5d\t%5d\t',ClusterMaxLocmm(:,i));
    fprintf(fid,'%10.2f\t%10d\n', ClusterMaxStat(i),ClusterSize(i));
end
% TODO: Write the data for all voxels to file also.

%% Do for the negative direction also
% find voxels above threshold
F = find(Img < -HeightThr);
% convert indices to voxel locations
DIM = Vin.dim;
[x y z] = ind2sub(DIM,F);
% HERE IS A USE OF AN SPM FUNCTION
% find clusters 
[A] = spm_clusters([x y z]');
% how many clusters
NumClust = max(A);
% how big are the clusters
ClusterSize = zeros(NumClust,1);
for i = 1:NumClust
    ClusterSize(i,1) = length(find(A == i));
end
% which clusters are big enough
SignClusters = find(ClusterSize >= ClusterThr);

SignVoxelsAndClusters = zeros(size(A));
ClusterMaxStat = zeros(length(SignClusters),1);
ClusterSize = zeros(length(SignClusters),1);
ClusterMaxLoc = [];
% what is the maximal value in the clusters and how big is each cluster
for i = 1:length(SignClusters)
    F = find(A == SignClusters(i));
    tempValues = zeros(size(F));
    ClusterSize(i) = length(tempValues);
    for j = 1:length(F)
        tempValues(j) = Img(x(F(j)), y(F(j)), z(F(j)));
    end
    ClusterMaxStat(i) = min(tempValues);
    tempMaxLoc = F(find(tempValues == ClusterMaxStat(i)));
    ClusterMaxLoc(i,:) = [x(tempMaxLoc) y(tempMaxLoc) z(tempMaxLoc)];
    SignVoxelsAndClusters(F) = 1;
end
% % if strmatch(upper(Flip),'Y')
% %     warning('Flipping results');
% %     for i = 1:size(ClusterMaxLoc,1)
% %         if ClusterMaxLoc(i,1)/floor(Vin.dim(1)/2) > 1
% %             ClusterMaxLoc(i,1)  = ClusterMaxLoc(i,1) - floor(Vin.dim(1)/2);
% %         elseif ClusterMaxLoc(i,1)/floor(Vin.dim(1)/2) < 1
% %             ClusterMaxLoc(i,1)  = ClusterMaxLoc(i,1) + floor(Vin.dim(1)/2);
% % 
% %         end
% %     end
% % end
ClusterSize
% find the locatio of the cluster maximum
ClusterMaxLocmm = Vin.mat*[ClusterMaxLoc ones(size(ClusterMaxLoc,1),1)]';
% ClusterMaxLocmm = ClusterMaxLocmm(1:3,:)';
% find the location of all significant voxels within significant clusters
xSign = x(find(SignVoxelsAndClusters));
ySign = y(find(SignVoxelsAndClusters));
zSign = z(find(SignVoxelsAndClusters));

XYZ = [xSign ySign zSign];
XYZmm = Vin.mat*[XYZ ones(length(XYZ),1)]';
XYZmm = XYZmm(1:3,:)';

% load the atlas maps
Vba = spm_vol(Pba);
Vaal = spm_vol(Paal);
Iaal = spm_read_vols(Vaal);
Iba = spm_read_vols(Vba);
% find the AAL and BA locations
%[AllAALList AllBAList] = subfnLocalFindAALandBA(XYZ, Iaal, Iba);
%[MaxAALList MaxBAList] = subfnLocalFindAALandBA(XYZ, Iaal, Iba);

[MaxAALList MaxBAList] = subfnLocalFindAALandBA(ClusterMaxLoc, Iaal, Iba);
for i = 1:size(ClusterMaxLocmm,2)
    fprintf(fid,'%20s\t%5s\t%5d\t',MaxAALList{i}(1:end-2), MaxAALList{i}(end),MaxBAList(i));
    fprintf(fid,'%5d\t%5d\t%5d\t',ClusterMaxLocmm(:,i));
    fprintf(fid,'%10.2f\t%10d\n', ClusterMaxStat(i),ClusterSize(i));
end

fclose(fid);
fprintf('Data saved to:\n\t%s\n',ClusterOutPutFileName);


    
function [AALList BAList] = subfnLocalFindAALandBA(XYZ, Iaal, Iba)
[aalCol1 aalCol2 aalCol3] = textread('/share/studies/CogRes/GroupAnalyses/ModMedCogRes/masks/aal.nii.txt','%d%s%d');
NVoxels = size(XYZ,1);
AALList = {};
BAList = zeros(NVoxels,1);
for i = 1:NVoxels
    CurrentValue = Iaal(XYZ(i,1), XYZ(i,2), XYZ(i,3));
    BAList(i,1) = Iba(XYZ(i,1), XYZ(i,2), XYZ(i,3));
    if CurrentValue
        AALList{i} = aalCol2{CurrentValue};
    else
        AALList{i} = '**empty**';
    end
end








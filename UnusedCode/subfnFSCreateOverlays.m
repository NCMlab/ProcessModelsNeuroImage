function subfnFSCreateOverlays(ResultsFolder, MaskList)
%[BasePath FileName] = fileparts(InputImage);
%BasePath = '/share/users/js2746_Jason/Studies/ModMedCogRes/Model58_XAgeGroup_MFreeSurfer_YFluid_WCogRes_VCogRes_COVnGMV_Sex_StudyID_Nboot5000_18-Mar-2013_11-11';
%P = fullfile(ResultsFolder,'BrainToParameters_JASON.nii');
P = '/Users/jason/Dropbox/SteffenerColumbia/Papers/MyPapers/Submitted/CogResMedMod/Figures/FREESURFEROVERLAYS/STANDARD_aparc+aseg.nii';

V = spm_vol(P);
I = spm_read_vols(V);
% Mask out a few regions
I(I==2) = 0;
I(I==41) = 0;
% P2 = InputImage;
% %P2 = fullfile(BasePath,'CondABMed1_pV0.00_sign0.0500_thickness.nii');
% %P2 = fullfile(BasePath,'Model2_FreeSurfer_t_thickness.nii');
% V2 = spm_vol(P2);
% I2 = spm_read_vols(V2);
% Load the mask image
% if nargin ~= 3
%     MaskList = unique(I);
% end
load('/Users/jason/Dropbox/SteffenerColumbia/Papers/MyPapers/Submitted/CogResMedMod/Figures/FREESURFEROVERLAYS/Mapping.mat');
Ibg = zeros(size(I));
BGList = [1:84];
 for i = 1:length(BGList)
     Ibg(find(I==Mapping(BGList(i)))) = Mapping(BGList(i));
 end
 I = Ibg;



% Create the mask
 Im = zeros(size(I));
 for i = 1:length(MaskList)
     Im(find(I==Mapping(MaskList(i)))) = 1;
 end
% thr = 1.96;
% % Take all statistical values and map them onto a range 
% posV = max(I2(find(I2>0)));
% negV = min(I2(find(I2<0)));
% 
% minRange = -3;
% maxRange = 15;
% %thr = 1.;
% posMap = hot(100);
% negMap = cool(100);
% 
[m n p] = size(I);
slices = [75 85 95 105 115 125 135 145 155];
%slices = [125];
Nslice = length(slices);
%%
% crop size
rect = [55 50 140 190];
Images = zeros(rect(4)+1,rect(3)+1,3,Nslice);
for i = 1:Nslice
    % Take the background image for this slice and rotate it
    BG = imrotate(squeeze(I(:,slices(i),:)),90);
    % dilate it so the cortical ribbon is a little wider which is better
    % for viewing
    BG =spm_dilate(BG);
    % invert the image
    BG = imcomplement(BG);
    % set all edge voxels to zero
    BG(find(BG<0))=0;
    BG = imcrop(BG,rect);
    % Load the mask if needed
    
    M = imrotate(squeeze(Im(:,slices(i),:)),90);
    M = imcrop(M,rect);
    Fm = M~=0;
    % set all three color channels to be the same and therefo, white
    outRGB1 = BG;
    outRGB2 = BG;
    outRGB3 = BG;
%     % load up the data slice
%     S1 = imrotate(squeeze(I2(:,slices(i),:)),90);
%    % S1 = spm_dilate(S1);
%     S1 = imcrop(S1,rect);
%     if ~isempty(Im)
%         S1(find(~Fm)) = 0;
%     end
%     
%     
%     FnegS1 = find(S1<-1.*thr);
%     FposS1 = find(S1>thr);
%     negS1 = S1(FnegS1)./min(S1(FnegS1));
%     posS1 = S1(FposS1)./max(S1(FposS1));
 
    % set the color of the data voxels
    outRGB1(Fm) = 1;
    outRGB2(Fm) = 0.8;
    %outRGB3(FnegS1) = negS1;
    % add this slice to the image montage
    Images(:,:,1,Nslice - i + 1) = outRGB1;
    Images(:,:,2,Nslice - i + 1) = outRGB2;
    Images(:,:,3,Nslice - i + 1) = outRGB3;
end
%colormap(hot(100))
f = figure;
montage(Images,'Size',[1 length(slices)])

[BasePath ResultsFolder] = fileparts(ResultsFolder);


set(f,'Name',ResultsFolder);

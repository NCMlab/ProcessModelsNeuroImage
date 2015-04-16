function WriteOutFSResultImages(Locations,Values,OutFile)

%BasePath = '/share/users/js2746_Jason/Studies/ModMedCogRes/Model58_XAgeGroup_MFreeSurfer_YFluid_WCogRes_VCogRes_COVnGMV_Sex_StudyID_Nboot5000_18-Mar-2013_11-11';
%P = fullfile(ResultsFolder,'BrainToParameters_JASON.nii');
if ismac
    BasePath = '/Users/jason/Dropbox/SteffenerColumbia/Scripts/ProcessModelsNeuroImage/FreeSurferFiles';
elseif isunix
    BasePath = '/home/jason/Dropbox/SteffenerColumbia/Scripts/ProcessModelsNeuroImage/FreeSurferFiles';
end

P = fullfile(BasePath, 'STANDARD_aparc+aseg.nii');

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
load(fullfile(BasePath,'Mapping.mat'));
Ibg = zeros(size(I));
BGList = [1:84];
 for i = 1:length(BGList)
     Ibg(find(I==Mapping(BGList(i)))) = Mapping(BGList(i));
 end
 I = Ibg;



% Create the mask
 Im = zeros(size(I));
 for i = 1:length(Locations)
     Im(find(I==Mapping(Locations(i)))) = Values(i);
 end
 V.fname = OutFile;
 V.dt = [16 0];
 spm_write_vol(V,Im);
 
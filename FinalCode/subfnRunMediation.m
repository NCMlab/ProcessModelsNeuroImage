function subfnRunMediation(Xname,X,Yname,Y,Mname,Mcontrast, subid,visitid, ResultsName, Nboot, ModelNum,InDir,OutDir)


NSubs = length(X);
SPMVer = 'spm8';
maskF = '/share/studies/iLS/GroupAnalyses/masks/mask_grey_p25.nii';
maskI = spm_read_vols(spm_vol(maskF));

VoxelsToInclude = find(maskI);
clear maskI;
Nvoxels = length(VoxelsToInclude);

IMGdata = zeros(NSubs,1,Nvoxels);

count = 1;
for i = 1:NSubs
        tempPath = fullfile(InDir,'Subjects',subid(i,:),visitid(i,end-4:end),'fmriStats',SPMVer,ResultsName,sprintf('con_%04d.img',Mcontrast));
        fprintf(1,'%s\n',tempPath);

        %tempPath = fullfile(BaseDir,'Subjects',subid{i},visitid{i}(end-4:end),'pASL',SPMVer,['swmeanCBF_0_mrfix_x1_pASL_' subid{i} '_' visitid{i}(end-4:end) '.img']);
        V = spm_vol(tempPath);
        temp = spm_read_vols(V);
        IMGdata(count,1,:) = temp(VoxelsToInclude);
        count = count + 1;
end

[m n p] = size(IMGdata);

AllData = {};
AllData.names = {};
AllData.names.Y = Yname;
AllData.names.X = Xname;
AllData.names.M = {Mname};
AllData.names.V = '';
AllData.names.W = '';
AllData.names.Q = '';
AllData.names.R = '';

AllData.Y = Y;
%AllData.M = reshape(IMGdata,m,1,n);
AllData.M = IMGdata;
AllData.X = X;
AllData.COV = [];
AllData.STRAT = X;
AllData.V = [];
AllData.W = [];
AllData.Q = [];
AllData.R = [];
AllData.Indices = VoxelsToInclude;
AllData.DIM = V.dim;
%Nboot = 5000;
%ModelNum = '4';
Thresholds = [0.05 0.01 0.005];
NJobSplit = 150;
[Nsub Nmed Nvoxels] = size(AllData.M);

% Create a parameter file
AnalysisParameters = {};
AnalysisParameters.BaseDir = OutDir;%'/share/studies/iLS/GroupAnalyses/100312_YngOld/Mediation';
AnalysisParameters.Nsub = Nsub;
AnalysisParameters.Nmed = Nmed;
AnalysisParameters.Nvoxels = Nvoxels;
AnalysisParameters.NJobSplit = NJobSplit;
AnalysisParameters.Nboot = Nboot;
AnalysisParameters.ModelNum = ModelNum;
AnalysisParameters.Thresholds = Thresholds;
AnalysisParameters.Indices = AllData.Indices;
AnalysisParameters.names = AllData.names;
AnalysisParameters.V = V;
AnalysisParameters.V.fname = '';
AnalysisParameters.V.descrip = '';
AnalysisParameters.Tag = ['X' Xname '_M' Mname '_Y' Yname];

subfnRunVoxelwiseProcess(AllData,AnalysisParameters);

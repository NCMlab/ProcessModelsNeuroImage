if ismac
    BaseDir = '/Users/jason/Dropbox/SteffenerColumbia/Scripts/ProcessModelsNeuroImage/TFCE';
elseif isunix
    BaseDir = '/home/jason/Dropbox/SteffenerColumbia/Scripts/ProcessModelsNeuroImage/TFCE';
end
setenv('FSLOUTPUTTYPE','NIFTI')

Vo = {};
Vo.fname = fullfile(BaseDir,'test1.nii');
Vo.dim = [20 20 20];
Vo.dt = [64 0];
Vo.mat = eye(4);
signal = zeros(Vo.dim);
signal(5:8,5:8,5:7) = -0.5;
signal(14:17,14:17,14:17) = -0.5;
signal(5:8,5:8,14:17) = -1;
Vs = Vo;
Vs.fname = fullfile(BaseDir,'signal.nii');
spm_write_vol(Vs,signal);

mask = zeros(Vo.dim);
mask(2:18,2:18,2:18) = 1;


Vm = Vo;
Vm.fname = fullfile(BaseDir,'mask.nii');
spm_write_vol(Vm,mask);

for i = 1:30
    Vo.fname = fullfile(BaseDir,sprintf('test_%04d.nii',i));
    I = signal + randn(Vo.dim);
    spm_write_vol(Vo,I);
end

t = tic;
I = randn(Vo.dim);



%setenv('FSLOUTPUTTYPE','NIFTI')
[a FSLMathsPath] = unix('! which fslmaths');
[a FSLStatsPath] = unix('! which fslstats');
[a FSLrandomisePath] = unix('! which randomise');
[a FSLMergePath] = unix('! which fslmerge');

inFile = Vo.fname;
outFile = fullfile(BaseDir,'test1out.nii');
Str = sprintf('! %s %s -mas %s -tfce 2 0.5 6 %s',FSLMathsPath(1:end-1),Vo.fname,Vm.fname,outFile)
unix(Str)
Str = sprintf('! %s %s -n -R',FSLStatsPath(1:end-1), outFile);
[a b] = unix(Str);
findUnder = find(b == ' ');
m = str2num(b(findUnder(1)+1:findUnder(2)-1))
toc(t)
Str = sprintf('! %s -t test4D test_0*',FSLMergePath(1:end-1));
unix(Str);

Str =sprintf(' ! %s -i test4D.nii.gz -o randomTest -d TFCE.mat -t TFCE.con --T2 -n 500 -m mask',FSLrandomisePath(1:end-1));
unix(Str)

%% run some stats 
% I have randomise results
% use the fslmaths  tfce function to compare to randomise
% Create the maximum statistic distribution

% Read data
Vi = spm_vol('test4D.nii');
I = spm_read_vols(Vi);
Nsub = 30;
Nperm = 500;
MaxStat = zeros(Nperm,1);
Indices = find(mask);
data = zeros(length(Indices),Nsub);
for i = 1:Nsub
    temp = I(:,:,:,i);
    data(:,i) = temp(Indices);
end
Io = zeros(Vo.dim);
outFile = fullfile(BaseDir,'test1out.nii');

for i = 1:Nperm
    pT = pinv([ones(Nsub,1) sign(randn(Nsub,1))]);
    temp = zeros(length(Indices),1);
    for j = 1:length(Indices)
        tempBeta = pT*data(j,:)';
        temp(j) = tempBeta(2);
    end
    % save the data to an image
    Io(Indices) = temp;
    Vo.fname = fullfile(BaseDir,'temp.nii');
    spm_write_vol(Vo,Io);
    % Run tfce and get the max value
    Str = sprintf('! %s %s -mas %s -tfce 2 0.5 6 %s',FSLMathsPath(1:end-1),Vo.fname,Vm.fname,outFile);
    unix(Str);
    Str = sprintf('! %s %s -n -R',FSLStatsPath(1:end-1),outFile);
    [a b] = unix(Str);
    findUnder = find(b == ' ');
    MaxStat(i) = str2num(b(findUnder(1)+1:findUnder(2)-1));
end
% Point estimate
pe = mean(I,4)+2;


% Run TFCE on this
Vo.fname = fullfile(BaseDir,'PointEstimate.nii');
outFile = fullfile(BaseDir,'PointEstimate_tfce.nii');
spm_write_vol(Vo,abs(pe));
% Run tfce and get the max value
Str = sprintf('! %s %s -mas %s -tfce 2 0.5 6 %s',FSLMathsPath(1:end-1),Vo.fname,Vm.fname,outFile);
unix(Str);
VpeTFCE = spm_vol(outFile);
IpeTFCE = spm_read_vols(VpeTFCE);
figure(2)
hist(IpeTFCE(Indices))

% Find the TFCE values of this and then use them as the point estimate
% For each voxel in the pe, compare it to the distribution of MaxStat values
pMap = zeros(Vo.dim);
for i = 1:length(Indices)
    pVal = sum(IpeTFCE(Indices(i)) > MaxStat)/Nperm;
    if pVal == 1
        pVal = 1/Nperm;
    elseif pVal == 0
        pVal = 1;
    end
    pMap(Indices(i)) = 1-pVal;
end
Vo.fname = fullfile(BaseDir,'pMap_tfce.nii');
spm_write_vol(Vo,pMap);
    


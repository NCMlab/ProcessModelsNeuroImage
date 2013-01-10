function BootData = subfnCompileClusterBootStrap(JobsDir)
% Findthe BootStrap results files
F = dir(fullfile(JobsDir,'BootStrapChunk*.mat'));
load(fullfile(JobsDir,F(1).name));
NChunks = length(F);
[NVoxels Nboot] = size(BootStrapResampleImages);
NTotalBoot = NChunks*Nboot;
BootData = zeros(NVoxels, NTotalBoot);

for i = 1:NChunks
    temp = load(fullfile(JobsDir,F(i).name));
    BootData(:,(i-1)*Nboot+1:i*Nboot) = temp.BootStrapResampleImages;
end
clear temp
% Clean up the files created
for i = 1:NChunks
    unix(sprintf('rm %s',fullfile(JobsDir,F(i).name)));
end
unix(['rm BOOTStrap_job.sh.*']);
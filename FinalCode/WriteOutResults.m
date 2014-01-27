function WriteOutResults(SelectedPath)
if nargin == 0
    SelectedPath = spm_select(1,'dir');
end
cd(SelectedPath)

if exist('AnalysisParameters.mat')
    load AnalysisParameters
else
    errordlg('this folder does not have the required AnalysticParameters.mat file');
end
% Check to see if the process is finished
F = dir('Results_*.mat');
if AnalysisParameters.NJobSplit == 1
    NvoxelsPerJob = ceil(AnalysisParameters.Nvoxels/AnalysisParameters.NJobSplit);
    % prealocate memory for AllParameters
    AllParameters = cell(AnalysisParameters.Nvoxels,1);
    % load up all the data

        clear Parameters
        load(F(1).name)
        AllParameters(:) = Parameters;
    

    [tVoxelIndices tImageVoxelIndices] = subfnWriteOutResults(AllParameters,AnalysisParameters,SelectedPath);

elseif AnalysisParameters.NJobSplit == length(F) 
    % the process is finished
    NvoxelsPerJob = ceil(AnalysisParameters.Nvoxels/AnalysisParameters.NJobSplit);
    % prealocate memory for AllParameters
    AllParameters = cell(AnalysisParameters.Nvoxels,1);
    % load up all the data
    for i = 1:length(F)-1
        clear Parameters
        load(F(i).name)
        AllParameters((i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob) = Parameters;
    end
    clear Parameters
    load(F(i+1).name)
    AllParameters(i*NvoxelsPerJob+1:end) = Parameters;

    [tVoxelIndices tImageVoxelIndices] = subfnWriteOutResults(AllParameters,AnalysisParameters,SelectedPath);
else
    errordlg('This process has not finished')
end

% for i = 1:length(AllParameters)
%     if ~isempty(AllParameters{i})
%         i
%         break
%     end
% end
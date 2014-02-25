if nargin == 0
    SelectedPath = spm_select(1,'dir');
end
cd(SelectedPath)
if exist('AnalysisParameters.mat')
    load AnalysisParameters
else
    errordlg('this folder does not have the required AnalysticParameters.mat file');
end
cd('Results')
% Check to see if the process is finished
F = dir('Results_*.mat');
NJobSplit = AnalysisParameters.NJobSplit;

NThr = length(AnalysisParameters.Thresholds);
index = 1;
Nvoxels = AnalysisParameters.Nvoxels;
OutData{index}.name = sprintf('CondM1M2_pV0');
OutData{index}.data = zeros(Nvoxels,1);
OutData{index}.field = [sprintf('M1M2{1}.p')];
OutData{index}.dataType = 16;
index = index + 1;
OutData{index}.name = sprintf('CondM1M2_pV1');
OutData{index}.data = zeros(Nvoxels,1);
OutData{index}.field = [sprintf('M1M2{2}.p')];
OutData{index}.dataType = 16;
index = index + 1;

if NJobSplit == length(F)
    % Load up the first chunk of data to be used to create the output data
    % structures
    NvoxelsPerJob = ceil(AnalysisParameters.Nvoxels/AnalysisParameters.NJobSplit);
    % Load just one chunk
    fprintf(1,'Loading data chunk %d of %d\n',1,NJobSplit);
    i = 1;
    clear Parameters
    load(F(i).name)
    IndicesForDataChunk = (i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob;
    Fields = fieldnames(Parameters{i}.M1M2{1}.BCaci);
    
    
    for i = 1:length(IndicesForDataChunk)
        for k = 1:2
            for j = 1:NThr
                if prod(getfield(Parameters{i}.M1M2{k}.BCaci,Fields{j}))>0
                    OutData{k}.data(IndicesForDataChunk(i)) = Parameters{i}.Thresholds(j);
                    break
                end
            end
        end
    end
    if NJobSplit > 1
        for i = 2:length(F)-1
            fprintf(1,'Loading data chunk %d of %d\n',i,NJobSplit);
            % Load the rest of the chunks one at a time
            clear Parameters
            load(F(i).name)
            IndicesForDataChunk = (i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob;
            for i = 1:length(IndicesForDataChunk)
                for k = 1:2
                    for j = 1:NThr
                        if prod(getfield(Parameters{i}.M1M2{k}.BCaci,Fields{j}))>0
                            OutData{k}.data(IndicesForDataChunk(i)) = Parameters{i}.Thresholds(j);
                            break
                        end
                    end
                end
            end
        end
    end

    
V = AnalysisParameters.V;
WIPsubfnWriteResultsToImages(OutData,AnalysisParameters.Indices,V,SelectedPath)

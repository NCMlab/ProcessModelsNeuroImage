function WIPWriteOutResults(SelectedPath)
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
    
    % prealocate memory for AllParameters
    OutData = subfnCreateOutputStructures(Parameters,AnalysisParameters,SelectedPath,0.05);
    % cyle over all the results and enter them into OutData
    OutData = subfnPutDataIntoOutputStructure(OutData,Parameters,AnalysisParameters,IndicesForDataChunk);
    if NJobSplit > 1
        for i = 2:length(F)-1
            fprintf(1,'Loading data chunk %d of %d\n',i,NJobSplit);
            % Load the rest of the chunks one at a time
            clear Parameters
            load(F(i).name)
            IndicesForDataChunk = (i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob;
            OutData = subfnPutDataIntoOutputStructure(OutData,Parameters,AnalysisParameters,IndicesForDataChunk);
        end
        fprintf(1,'Loading data chunk %d of %d\n',i+1,NJobSplit);
        clear Parameters
        load(F(i+1).name)
        IndicesForDataChunk = i*NvoxelsPerJob+1:AnalysisParameters.Nvoxels;
        OutData = subfnPutDataIntoOutputStructure(OutData,Parameters,AnalysisParameters,IndicesForDataChunk);


        
    else
        errordlg('This process has not finished')
    end
end    



 
%%
% Write out the images
V = AnalysisParameters.V;
WIPsubfnWriteResultsToImages(OutData,AnalysisParameters.Indices,V,SelectedPath)


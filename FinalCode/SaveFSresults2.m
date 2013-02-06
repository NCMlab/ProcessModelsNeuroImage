function SaveFSresults2
%BasePath = '/Users/jason/Documents/MyData/ModMedCogRes';
%BasePath = 'C:\Users\steffener\Dropbox\ModMedCogRes';
BasePath = 'w:/js2746_Jason/Scripts/ProcessModelsNeuroImage/FreeSurferFiles';
BasePath = '/share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FreeSurferFiles';

load(fullfile(BasePath,'FSheader.mat'));
SelectedPath = spm_select(1,'dir');
cd(SelectedPath)
if exist('AnalysisParameters.mat')
    load AnalysisParameters
else
    errordlg('this folder does not have the required AnalysticParameters.mat file');
end

MeasureOfInterest = 'volume';
% Check to see if the process is finished
F = dir('Results_*.mat');
if AnalysisParameters.NJobSplit == length(F) 
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
    load(F(1).name)
    AllParameters = Parameters;
    
     OutData = subfnCreateOutDataStructureForModels(AllParameters,AnalysisParameters);
     OutData = subfnAddConditionalEffectsToOutPut(OutData,AnalysisParameters);
     OutData = subfnPutDataIntoOutputStructure(OutData,AllParameters,AnalysisParameters);
     subfnWriteFSResultsToNIFTI(OutData,SelectedPath,BasePath,MeasureOfInterest) 
      
else
    errordlg('This process has not finished')
end

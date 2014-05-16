function WriteOutFreeSurferResults
% This function writes out results to CSV files instead of images.
% As of rigth now it does not write out conditional effects
BasePath = '/Users/jason/Documents/MyData/ModMedCogRes';
BasePath = '/share/users/js2746_Jason/Studies/ModMedCogRes/DataFiles';
load(fullfile(BasePath,'FSheader.mat'));

SelectedPath = spm_select(1,'dir');
cd(SelectedPath)

if exist('AnalysisParameters.mat')
    load AnalysisParameters
else
    errordlg('this folder does not have the required AnalysticParameters.mat file');
end
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
    load(F(end).name);
    AllParameters((length(F)-1)*NvoxelsPerJob+1:end) = Parameters;
    

     [OutData index] = subfnCreateOutDataStructureForModels(AllParameters,AnalysisParameters);
        OutData = subfnPutDataIntoOutputStructure(OutData,AllParameters,AnalysisParameters)
else
    errordlg('This process has not finished')
end


% for i = 1:length(AllParameters)
%     if ~isempty(AllParameters{i})
%         i
%         break
%     end
% end


%%
HeaderFlag = 1;
fid = 1;
fid = fopen('FSResults.csv','w');
fprintf(fid,'%30s,','RegionName');
for j = 1:length(OutData)
    fprintf(fid,'%30s,',OutData{j}.name);
end
fprintf(fid,'\n');
for i = 1:NvoxelsPerJob
    fprintf(fid,'%-30s,',Header{i});
    for j = 1:length(OutData)
        fprintf(fid,'%30.4f,',OutData{j}.data(i));
    end
    fprintf(fid,'\n');
end
fclose(fid);
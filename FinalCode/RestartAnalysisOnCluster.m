function RestartAnalysisOnCluster
% If a model being estimated on the cluster was stopped and not finished.
% This program resubmits the unfinished jobs to the cluster and uses the
% existing folders.
%
SelectedPath = spm_select(1,'dir');
cd(SelectedPath)

OutFolder = SelectedPath;
JobFolder = fullfile(OutFolder,'jobs');

if exist('AnalysisParameters.mat')
    load AnalysisParameters
else
    errordlg('this folder does not have the required AnalysticParameters.mat file');
end
% Check to see if the process is finished
F = dir('Results_*.mat');
if ~(AnalysisParameters.NJobSplit == length(F))
    % find the unfinished jobs
    FinishedData = zeros(AnalysisParameters.NJobSplit,1);
    for i = 1:AnalysisParameters.NJobSplit
        if ~isempty(dir(sprintf('Results_%04d.mat',i)))
            FinishedData(i) = 1;
        end
    end
    for i = 1:AnalysisParameters.NJobSplit
        if ~FinishedData(i)
            jobPath = fullfile(OutFolder,sprintf('job_%04d.sh',i));
            %    Str = ['! qsub  ' jobPath];
            Str = ['! qsub -q short.q -p -10 -e ' JobFolder ' -o ' JobFolder ' -l mem_free=500M ' jobPath];
            unix(Str);
        end
    end
end
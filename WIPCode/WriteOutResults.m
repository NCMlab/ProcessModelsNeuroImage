function WriteOutResults(ResultsPath)
% Take the results from an analysis and write out results as NIFTI images.
if nargin == 0
    ResultsPath = spm_select(1,'dir','Select analysis directory');
    % Check to make sure this seems correct
    % To Do
end

% Load the data/parameters used in this analysis
load(fullfile(ResultsPath,'data','ModelInfo'))
ModelInfo.ResultsPath = ResultsPath;

% Remove the data from the structure to preserve memory
ModelInfo.data = [];

% Check to see if there are Permute results
if ~isempty(dir(fullfile(ResultsPath,'Results','Permute*.mat')))
    ModelType = 'permutation';
elseif ~isempty(dir(fullfile(ResultsPath,'Results','Bootstrap*.mat')))
    ModelType = 'bootstrap';
else
    errordlg('Unknown results');
end

% Find the number of results files
switch ModelType
    case 'permutation'
        % locate the results files
        F = dir(fullfile(ResultsPath,'Results','Permute*.mat'));
        NFiles = length(F);
        % load a single results fle to determine the size of the paths
        load(fullfile(ResultsPath,'Results',F(1).name))
        [m n o p] = size(MaxPaths);
        % Check to make sure all the files are there
        if ~(ModelInfo.Nperm == p*NFiles)
            errordlg('There are not enough results files based onthe specified parameters.');
        end
        % Create the structures to hold the permutation results for the
        % path values and the standardized parameter estimates
        MaxPermPaths = zeros(m,n,o,p*NFiles);
        MinPermPaths = zeros(m,n,o,p*NFiles);
        [mB nB pB] = size(MaxB);
        MaxPermB = zeros(mB,nB,pB*NFiles);
        MinPermB = zeros(mB,nB,pB*NFiles);
        
        % load the data and put the permutation results in these structures
        for i = 1:NFiles
            load(fullfile(ResultsPath,'Results',F(i).name))
            MaxPermPaths(:,:,:,(i-1)*p+1:i*p) = MaxPaths;
            MinPermPaths(:,:,:,(i-1)*p+1:i*p) = MinPaths;
            MaxPermB(:,:,(i-1)*p+1:i*p) = MaxB;
            MinPermB(:,:,(i-1)*p+1:i*p) = MinB;
        end
        % Load up the point estimate results
        F = dir(fullfile(ResultsPath,'Results','PointEstimate*.mat'));
        load(fullfile(ResultsPath,'Results',F(1).name))
        
        % I am not sure what this is supposed to do
                PointEstimate = zeros(m,n,o,length(Parameters));
                for i = 1:length(Parameters)
                    for kk = 1:o
                        PointEstimate(:,:,kk,i) = Parameters{i}.Paths{kk};
                    end
                end
    case 'bootstrap'
end

%% WRITE OUT ALL IMAGES
WriteOutParameterMaps('beta',Parameters,ModelInfo)
WriteOutParameterMaps('B',Parameters,ModelInfo)
WriteOutParameterMaps('t',Parameters,ModelInfo)


%% PATH IMAGES

%load(fullfile(fileparts(pwd),'data','AllData'));

for j = 1:length(ModelInfo.Thresholds)
    c = floor(ModelInfo.Thresholds(j)*ModelInfo.Nperm);
    if c == 0
        c = 1;
    end
    
    for kk = 1:o % cycle over the number of paths
        for i = 1:m
            sMax(:,i) = sort(squeeze(MaxPermPaths(i,:,kk,:)),'descend');
            sMin(:,i) = sort(squeeze(MinPermPaths(i,:,kk,:)));
            Mx(i) = sMax(c,i);
            Mn(i) = sMin(c,i);
            % write unthresholed path estimate
            Vo = ModelInfo.DataHeader;
            Vo.fname = (fullfile(ModelInfo.ResultsPath,sprintf('Path%d_level%d.nii',kk,i)));
            I = zeros(ModelInfo.DataHeader.dim);
            temp = squeeze(PointEstimate(i,1,kk,:));
             
            
            I(ModelInfo.Indices) = temp;
            spm_write_vol(Vo,I);
            
            temp(find((PointEstimate(i,1,kk,:) < Mx(i))&(PointEstimate(i,1,kk,:) >0))) = 0;
            temp(find((PointEstimate(i,1,kk,:) > Mn(i))&(PointEstimate(i,1,kk,:) <0))) = 0;
            
            
            
            Vo.fname=(fullfile(ModelInfo.ResultsPath,sprintf('Path%d_level%d_%0.4f.nii',kk,i,ModelInfo.Thresholds(j))));
            I = zeros(Vo.dim);
            I(ModelInfo.Indices) = temp;
            spm_write_vol(Vo,I);
            
        end
    end
end


function WriteOutResults(ResultsFolder)
% Take the results from an analysis and write out results as NIFTI images.
if nargin == 0
    ResultsFolder = spm_select(1,'dir','Select analysis directory');
end

% Load the data/parameters used in this analysis
load(fullfile(ResultsFolder,'data','ModelInfo'))
ModelInfo.ResultsPath = ResultsFolder;

% Remove the data from the structure to preserve memory
ModelInfo.data = [];

% Check to see if there are Permute results
if ~isempty(dir(fullfile(ResultsFolder,'Results','Permute*.mat')))
    ModelType = 'permutation';
elseif ~isempty(dir(fullfile(ResultsFolder,'Results','BootStrap*.mat')))
    ModelType = 'bootstrap';
else
    errordlg('Unknown results');
end

% Find the number of results files
switch ModelType
    case 'permutation'
        % locate the results files
        F = dir(fullfile(ResultsFolder,'Results','Permute*.mat'));
        NFiles = length(F);
        
        % Check to see if the analyses are completed
        if ~(NFiles == ModelInfo.NJobSplit)
            errordlg('This analysis is not complete.');
        end
        
        % load a single results fle to determine the size of the paths
        load(fullfile(ResultsFolder,'Results',F(1).name))
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
            load(fullfile(ResultsFolder,'Results',F(i).name))
            MaxPermPaths(:,:,:,(i-1)*p+1:i*p) = MaxPaths;
            MinPermPaths(:,:,:,(i-1)*p+1:i*p) = MinPaths;
            MaxPermB(:,:,(i-1)*p+1:i*p) = MaxB;
            MinPermB(:,:,(i-1)*p+1:i*p) = MinB;
        end
        % Load up the point estimate results
        F = dir(fullfile(ResultsFolder,'Results','PointEstimate*.mat'));
        load(fullfile(ResultsFolder,'Results',F(1).name))
        
        % Reshape the path estimates for the next step of writing them out
        % as images
        PointEstimate = zeros(m,n,o,length(Parameters));
        for i = 1:length(Parameters)
            for kk = 1:o
                PointEstimate(:,:,kk,i) = Parameters{i}.Paths{kk};
            end
        end
        
        WriteOutPermutationPaths(ModelInfo,MaxPermPaths,MinPermPaths,PointEstimate,o,m)
    case 'bootstrap'
        % There is a problem here fusing the split results back together
        
        % locate the results files
        F = dir(fullfile(ResultsFolder,'Results','BootStrap*.mat'));
        NFiles = length(F);
        
        % Check to see if the analyses are completed
        if ~(NFiles == ModelInfo.NJobSplit)
            errordlg('This analysis is not complete.');
        end
        
        % load a single results fle to determine the size of the paths
        load(fullfile(ResultsFolder,'Results',F(1).name))
        [n m] = size(Results{1}.Paths{1});
        PointEstimate = zeros(m,n,ModelInfo.Nvoxels);
        % cycle over number of voxels
        Parameters = cell(ModelInfo.Nvoxels,1);
        for k = 1:NFiles
            load(fullfile(ResultsFolder,'Results',F(k).name))
            for i = 1:length(Results)
                Parameters{(k-1)*ModelInfo.NJobSplit + i} = Results{i};
                PointEstimate(:,:,(k-1)*ModelInfo.NJobSplit + i) = Results{i}.Paths{:};
                % cycle over thresholds
                for j = 1:length(ModelInfo.Thresholds)
                    BCaCIUpper(:,:,j,(k-1)*ModelInfo.NJobSplit + i) = Results{i}.BCaCI.Paths(:,:,j,1,1);
                    BCaCILower(:,:,j,(k-1)*ModelInfo.NJobSplit + i) = Results{i}.BCaCI.Paths(:,:,j,1,2);
                end
            end
        end
        WriteOutBootstrapPaths(ModelInfo,PointEstimate,BCaCIUpper,BCaCILower,m,n)
end






%% WRITE OUT ALL IMAGES from the regression models
WriteOutParameterMaps('beta',Parameters,ModelInfo)
WriteOutParameterMaps('B',Parameters,ModelInfo)
WriteOutParameterMaps('t',Parameters,ModelInfo)

%%
function WriteOutBootstrapPaths(ModelInfo,PointEstimate,BCaCIUpper,BCaCILower,m,n)
%% WRITE OUT THE PATH IMAGES FOR THE BOOTSTRAP TEST

for i = 1:m % cycle over path number
    for j = 1:n % cycle over probed level in the path
        
        % write unthresholed path estimate
        Vo = ModelInfo.DataHeader;
        Vo.fname = (fullfile(ModelInfo.ResultsPath,sprintf('Path%d_level%d.nii',i,j)));
        % Prepare the data matrix
        I = zeros(ModelInfo.DataHeader.dim);
        % Extract the path data values
        temp = squeeze(PointEstimate(i,j,:));
        I(ModelInfo.Indices) = temp;
        % Write the images
        spm_write_vol(Vo,I);
    end
end
% Write out the thresholded path images
for k = 1:length(ModelInfo.Thresholds)
    for i = 1:m % cycle over path number
        for j = 1:n % cycle over probed level in the path
            
            % Find those locations where the product ofthe confidence
            % intervals is greater then zero. These are locatons where the
            % CI do not include zero.
            temp = squeeze(BCaCILower(j,i,k,:).*BCaCIUpper(j,i,k,:)) > 0;
            Vo = ModelInfo.DataHeader;
            Vo.fname = (fullfile(ModelInfo.ResultsPath,sprintf('Path%d_level%d_%0.3f.nii',i,j,ModelInfo.Thresholds(k))));
            % Prepare the data matrix
            I = zeros(ModelInfo.DataHeader.dim);

            I(ModelInfo.Indices) = temp;
            % Write the images
            spm_write_vol(Vo,I);
        end
    end
end



function WriteOutPermutationPaths(ModelInfo,MaxPermPaths,MinPermPaths,PointEstimate,o,m)
%% WRITE OUT THE PATH IMAGES FOR THE PERMUTATION TEST
% cycle over the thresholds requested
for j = 1:length(ModelInfo.Thresholds)
    % Find the number in a sorted list of permutations that corresponds to
    % the threshold.
    c = floor(ModelInfo.Thresholds(j)*ModelInfo.Nperm);
    % If the value exceeds the precision of the number of permutaions then
    % set it to zero.
    % e.g. a threshold of 0.00001 is not possible with 100 permutations.
    % The most precise threshold is: 1/100 = 0.01
    if c == 0
        % RESET the threshold used to the most precise
        ModelInfo.Thresholds(j) = 1/ModelInfo.Nperm;
        c = 1;
    end
    
    for kk = 1:o % cycle over the number of paths
        % cycle over the probed value??
        for i = 1:m
            % sort the max and min permutation results
            sMax(:,i) = sort(squeeze(MaxPermPaths(i,:,kk,:)),'descend');
            sMin(:,i) = sort(squeeze(MinPermPaths(i,:,kk,:)));
            % find the permutation value based on the sorted values
            Mx(i) = sMax(c,i);
            Mn(i) = sMin(c,i);
            
            % write unthresholed path estimate
            Vo = ModelInfo.DataHeader;
            Vo.fname = (fullfile(ModelInfo.ResultsPath,sprintf('Path%d_level%d.nii',kk,i)));
            % Prepare the data matrix
            I = zeros(ModelInfo.DataHeader.dim);
            % Extract the path data values
            temp = squeeze(PointEstimate(i,1,kk,:));
            I(ModelInfo.Indices) = temp;
            % Write the images
            spm_write_vol(Vo,I);
            
            % Find the locations that exceed the two-tailed threshold
% >>>> I THINK THERE IS A BUG HERE WITH THE SECOND INDEX BEING SET TO 1 <<<
            temp(find((PointEstimate(i,1,kk,:) < Mx(i))&(PointEstimate(i,1,kk,:) >0))) = 0;
            temp(find((PointEstimate(i,1,kk,:) > Mn(i))&(PointEstimate(i,1,kk,:) <0))) = 0;
            
            % Create the thresholded path output file name.
            Vo.fname=(fullfile(ModelInfo.ResultsPath,sprintf('Path%d_level%d_%0.4f.nii',kk,i,ModelInfo.Thresholds(j))));
            I = zeros(Vo.dim);
            I(ModelInfo.Indices) = temp;
            % Write the image.
            spm_write_vol(Vo,I);
        end
    end
end


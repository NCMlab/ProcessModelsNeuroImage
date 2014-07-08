function WriteOutResults(ResultsFolder)
% Take the results from an analysis and write out results as NIFTI images.
if nargin == 0
    ResultsFolder = spm_select(1,'dir','Select analysis directory');
end

fprintf(1,'%s\n',ResultsFolder);

% Check to see if there are Permute results
if ~isempty(dir(fullfile(ResultsFolder,'Results','Permute*.mat')))
    ModelType = 'permutation';
elseif ~isempty(dir(fullfile(ResultsFolder,'Results','BootStrap*.mat')))
    ModelType = 'bootstrap';
else
    errordlg('Unknown results');
end
load(fullfile(ResultsFolder,'data','ModelInfo'))
ModelInfo.ResultsPath = ResultsFolder;

% Remove the data from the structure to preserve memory
ModelInfo.data = [];



% Find the number of results files
switch ModelType
    case 'permutation'
        % Load the data/parameters used in this analysis
       
        
        % locate the results files
        F = dir(fullfile(ResultsFolder,'Results','Permute*.mat'));
        NFiles = length(F);
        
        % Check to see if the analyses are completed
        if ~(NFiles == ModelInfo.NJobSplit)
            errordlg('This analysis is not complete.');           
        end
        
        % load a single results fle to determine the size of the paths
        load(fullfile(ResultsFolder,'Results',F(1).name))
        % m: steps in this moderated path for dimension 1
        % n: steps in this moderated path for dimension 2
        % o: number of paths
        % p: number of permutaions for this chunk of results
        [m n o p] = size(MaxPaths);
        % Check to make sure all the files are there
        if ~(ModelInfo.Nperm == p*NFiles)
            errordlg('There are not enough results files based onthe specified parameters.');
        end
        % Create the structures to hold the permutation results for the
        % path values and the standardized parameter estimates
        MaxPermPaths = zeros(m,n,o,p*NFiles);
        MinPermPaths = zeros(m,n,o,p*NFiles);
        [mB nB pB] = size(MaxBeta);
        MaxPermB = zeros(mB,nB,pB*NFiles);
        MinPermB = zeros(mB,nB,pB*NFiles);
        
        % load the data and put the permutation results in these structures
        for i = 1:NFiles
            load(fullfile(ResultsFolder,'Results',F(i).name))
            MaxPermPaths(:,:,:,(i-1)*p+1:i*p) = MaxPaths;
            MinPermPaths(:,:,:,(i-1)*p+1:i*p) = MinPaths;
            MaxPermB(:,:,(i-1)*p+1:i*p) = MaxBeta;
            MinPermB(:,:,(i-1)*p+1:i*p) = MinBeta;
        end
        % Load up the point estimate results
        F = dir(fullfile(ResultsFolder,'Results','PointEstimate*.mat'));
        load(fullfile(ResultsFolder,'Results',F(1).name))

        % Reshape the path estimates for the next step of writing them out
        % as images
        PointEstimatePath = zeros(m,n,o,length(Parameters));
        PointEstimateB = zeros(mB, nB, length(Parameters));
        for i = 1:length(Parameters)
            for kk = 1:o
                PointEstimatePath(:,:,kk,i) = Parameters{i}.Paths{kk};
            end
            PointEstimateB(:,:,i) = Parameters{i}.B;
        end
        
         % Display point estimate distribution
%         figure      
%         hist(PointEstimatePath(2,:,:,:),50)
%         
        
        AllParameters = Parameters;
        clear Parameters;
        WriteOutPermutationPaths(ModelInfo,MaxPermPaths,MinPermPaths,PointEstimatePath,o,m)
        WriteOutPermutationB(ModelInfo,MaxPermB,MinPermB,PointEstimateB)
        
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
        if exist('Results','var')
            Parameters = Results;
        end
        [n m] = size(Parameters{1}.Paths{1});
        % Prespecify the data structures
        PointEstimate = zeros(m,n,ModelInfo.Nvoxels);
        % Create a structure to contain the parameters from all analysis chunks
        AllParameters = cell(ModelInfo.Nvoxels,1);
        % m: path number
        % n: probed level in the path
        BCaCIUpper = zeros(n,m,length(ModelInfo.Thresholds),ModelInfo.Nvoxels);
        BCaCILower = zeros(n,m,length(ModelInfo.Thresholds),ModelInfo.Nvoxels);
        
        NvoxelsPerJob = ceil(ModelInfo.Nvoxels/ModelInfo.NJobSplit);
        
        Indices = [];
        % Cycle over each file
        for k = 1:NFiles
            % get the indices
              % load each results file
              load(fullfile(ResultsFolder,'Results',F(k).name))
              if exist('Results','var')
                  Parameters = Results;
              end

            %ModelInfo.Indices = [ModelInfo.Indices Parameters
            % cycle over the voxels in the results file
            for i = 1:length(Parameters)
                % what is the overall index of voxels as described by this
                % chunk of results
                
                Index = (k-1)*NvoxelsPerJob + i;
                AllParameters{Index} = Parameters{i};
                PointEstimate(:,:,Index) = Parameters{i}.Paths{:};
                % cycle over thresholds
                for j = 1:length(ModelInfo.Thresholds)
                    BCaCIUpper(:,:,j,Index) = Parameters{i}.BCaCI.Paths(:,:,1,:,j);
                    BCaCILower(:,:,j,Index) = Parameters{i}.BCaCI.Paths(:,:,2,:,j);
                end
            end
        end
        WriteOutBootstrapPaths(ModelInfo,PointEstimate,BCaCIUpper,BCaCILower,m,n)
        WriteOutParameterMaps('BCaCI.p',AllParameters,ModelInfo)
        

        % Write out the FDR thresholded p maps also

        WriteOutParameterMaps('BCaCI.p',AllParameters,ModelInfo,1)
        WriteOutParameterMaps('BCaCI.Z',AllParameters,ModelInfo)
         WriteOutSingleMap('BCaCI.PathsZ',AllParameters,ModelInfo)
         WriteOutSingleMap('BCaCI.PathsP',AllParameters,ModelInfo)
         WriteOutSingleMap('BCaCI.PathsP',AllParameters,ModelInfo,1)
end

%% WRITE OUT ALL IMAGES from the regression models
WriteOutParameterMaps('beta',AllParameters,ModelInfo)
WriteOutParameterMaps('B',AllParameters,ModelInfo)
WriteOutParameterMaps('t',AllParameters,ModelInfo)


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
% This allows the images to be written out even if the processing is not
% complete;
Nperm = size(MaxPermPaths,4);


for j = 1:length(ModelInfo.Thresholds)
    % Find the number in a sorted list of permutations that corresponds to
    % the threshold.
    c = floor(ModelInfo.Thresholds(j)*Nperm);
    % If the value exceeds the precision of the number of permutaions then
    % set it to zero.
    % e.g. a threshold of 0.00001 is not possible with 100 permutations.
    % The most precise threshold is: 1/100 = 0.01
    if c == 0
        % RESET the threshold used to the most precise
        ModelInfo.Thresholds(j) = 1/Nperm;
        c = 1;
    end
    
    for kk = 1:o % cycle over the number of paths
        
        % cycle over the probed value for dimension 1
        % The code does not handle multidimensional interactions yet.
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
            
            temp(find((PointEstimate(i,1,kk,:) < Mx(i)) & (PointEstimate(i,1,kk,:) >0))) = 0;
            temp(find((PointEstimate(i,1,kk,:) > Mn(i)) & (PointEstimate(i,1,kk,:) <0))) = 0;
            
            % Create the thresholded path output file name.
            Vo.fname=(fullfile(ModelInfo.ResultsPath,sprintf('Path%d_level%d_%0.4f.nii',kk,i,ModelInfo.Thresholds(j))));
            I = zeros(Vo.dim);
            I(ModelInfo.Indices) = temp;
            % Write the image.
            spm_write_vol(Vo,I);
        end
    end
end

function WriteOutPermutationB(ModelInfo,MaxB,MinB,PointEstimate)

% cycle over the thresholds requested
for j = 1:length(ModelInfo.Thresholds)
    Tag = sprintf('%s_%0.4f','Bperm',ModelInfo.Thresholds(j));
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
    % cycle over columns
    for i = 1:ModelInfo.Nvar
        if sum(ModelInfo.Direct(:,i))
            % CONSTANT TERMS ARE ROW 1
            
            % cycle over the rows in the model
            for j = 1:ModelInfo.Nvar
                if ModelInfo.Direct(j,i)
                    % create the filename describing the dependent and
                    % independent effects
                    FileName = sprintf('Model%d_DEP%s_IND%s_%s.nii',i,ModelInfo.Names{i},ModelInfo.Names{j},Tag);
                    % sort the max and min permutation results
                    sMax = sort(squeeze(MaxB(j+1,i,:)),'descend');
                    sMin = sort(squeeze(MinB(j+1,i,:)));
                    % find the permutation value based on the sorted values
                    Mx = sMax(c);
                    Mn = sMin(c);
                    
                    temp = squeeze(PointEstimate(j+1,i,:));
                    temp(find((PointEstimate(j+1,i,:) < Mx)&(PointEstimate(j+1,i,:) >0))) = 0;
                    temp(find((PointEstimate(j+1,i,:) > Mn)&(PointEstimate(j+1,i,:) <0))) = 0;
                    
                    I = zeros(ModelInfo.DataHeader.dim);
                    I(ModelInfo.Indices) = temp;
                    % Create the header for this image
                    Vo = ModelInfo.DataHeader;
                    Vo.fname = fullfile(ModelInfo.ResultsPath,FileName);
                    spm_write_vol(Vo,I);
                end
            end
            % Interaction terms
            if sum(ModelInfo.Inter(:,i))
                InterVar = find(ModelInfo.Inter(:,i));
                FileName = sprintf('Model%d_DEP%s_IND',i,ModelInfo.Names{i});
                for j = 1:length(InterVar)
                    FileName = sprintf('%s%sX',FileName,ModelInfo.Names{InterVar(j)});
                end
                FileName = sprintf('%s_%s.nii',FileName(1:end-1),Tag);

                % sort the max and min permutation results
                sMax = sort(squeeze(MaxB(ModelInfo.Nvar+2,i,:)),'descend');
                sMin = sort(squeeze(MinB(ModelInfo.Nvar+2,i,:)));
                % find the permutation value based on the sorted values
                Mx = sMax(c);
                Mn = sMin(c);
                
                temp = squeeze(PointEstimate(ModelInfo.Nvar+2,i,:));
                temp(find((PointEstimate(ModelInfo.Nvar+2,i,:) < Mx)&(PointEstimate(ModelInfo.Nvar+2,i,:) >0))) = 0;
                temp(find((PointEstimate(ModelInfo.Nvar+2,i,:) > Mn)&(PointEstimate(ModelInfo.Nvar+2,i,:) <0))) = 0;
                
                
                
                I = zeros(ModelInfo.DataHeader.dim);
                I(ModelInfo.Indices) = squeeze(temp);
                % Create the header for this image
                Vo = ModelInfo.DataHeader;
                Vo.fname = fullfile(ModelInfo.ResultsPath,FileName);
                spm_write_vol(Vo,I);
            end
        end
    end
end

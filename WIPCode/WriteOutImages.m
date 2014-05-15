function WriteOutResults(ResultsPath)
% Take the results from an analysis and write out results as NIFTI images.
if nargin == 0
    ResultsPath = spm_select(1,'dir','Select analysis directory');
    % Check to make sure this seems correct
    % To Do
end

% Load the data/parameters used in this analysis
load(fullfile(fileparts(ResultsPath),'data','ModelInfo'))

% Remove the data from the structure to preserve memory
AnalysisParameters = ModelInfo;
AnalysisParameters.data = [];

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
        load(F(1).name)
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
            load(F(i).name)
            MaxPermPaths(:,:,:,(i-1)*p+1:i*p) = MaxPaths;
            MinPermPaths(:,:,:,(i-1)*p+1:i*p) = MinPaths;
            MaxPermB(:,:,(i-1)*p+1:i*p) = MaxB;
            MinPermB(:,:,(i-1)*p+1:i*p) = MinB;
        end
        % Load up the point estimate results
        F = dir(fullfile(ResultsPath,'Results','PointEstimate*.mat'));
        load(F(1).name)
        PointEstimate = zeros(m,n,o,length(Parameters));
        for i = 1:length(Parameters)
            for kk = 1:o
                PointEstimate(:,:,kk,i) = Parameters{i}.Paths{kk};
            end
        end
    case 'bootstrap'
end

%% WRITE OUT ALL IMAGES 
WriteOutParameterMaps('beta')
WriteOutParameterMaps('B')
WriteOutParameterMaps('t')
%% BETA IMAGES PERMUTATION STATISTIC
Nvoxels = length(Parameters);
Tag = 'beta';
temp = zeros([size(Parameters{1}.beta) Nvoxels]);
for i =1:Nvoxels
    temp(:,:,i) = getfield(Parameters{i},Tag);
end
Nvar = AnalysisParameters.Nvar;
% 
% 
% 
% for k = 1%:length(ModelInfo.Thresholds)
%     c = floor(ModelInfo.Thresholds(k)*Nperm);
%     for i = 1:Nvar % Columns
%         % Is there something in this column?
%         if sum(ModelInfo.Direct(:,i))
%             for j = 1:Nvar % Rows
%                 if ModelInfo.Direct(j,i)
%                     tempPERM = squeeze(MaxPermB(j+1,i,:));
%                     tempMaxB = sort(tempPERM,'descend');
%                     tempMaxB = tempMaxB(c);
%                     figure
%                     clf
%                     subplot(132)
%                     hist(tempPERM,40)
%                     h = line([tempMaxB tempMaxB],[0 100]);
%                     set(h,'Color','r')
%                     
%                     tempPERM = squeeze(MinPermB(j+1,i,:));
%                     tempMinB = sort(tempPERM,'ascend');
%                     tempMinB = tempMinB(c);
%                     
%                     subplot(131)
%                     hist(tempPERM,40)
%                     h = line([tempMinB tempMinB],[0 100]);
%                     set(h,'Color','r')
%                     
%                     tempBETA = squeeze(temp(j+1,i,:));
%                     subplot(133)
%                     hist(tempBETA,40)
%                     h = line([tempMinB tempMinB],[0 1000]);
%                     set(h,'Color','r')
%                     h = line([tempMaxB tempMaxB],[0 1000]);
%                     set(h,'Color','r')
%                     
%                     tempBETA(find((tempBETA < tempMaxB) & (tempBETA > 0))) = 0;
%                     tempBETA(find((tempBETA > tempMinB) & (tempBETA < 0))) = 0;
%                     
%                     FileName = sprintf('Model%d_DEP%s_IND%s_%s_%0.4f.nii',i,ModelInfo.Names{i},ModelInfo.Names{j},Tag,ModelInfo.Thresholds(k));
%                     
%                     subplot(132)
%                     h = title(FileName);
%                     set(h,'Interpreter','none');
%                     
%                     I = zeros(V.dim);
%                     I(ModelInfo.Indices) = tempBETA;
%                     Vo = V;
%                     Vo.fname = fullfile(fileparts(pwd),FileName);
%                     spm_write_vol(Vo,I);
%                     
%                     
%                     % Interaction terms
%                     if sum(ModelInfo.Inter(:,i))
%                         InterVar = find(ModelInfo.Inter(:,i));
%                         FileName = sprintf('Model%d_DEP%s_IND',i,ModelInfo.Names{i});
%                         for j = 1:length(InterVar)
%                             FileName = sprintf('%s%sX',FileName,ModelInfo.Names{InterVar(j)});
%                         end
%                         FileName = sprintf('%s_%s.nii',FileName(1:end-1),Tag);
%                         
%                         tempPERM = squeeze(MaxPermB(j+1,i,:));
%                         tempMaxB = sort(tempPERM,'descend');
%                         tempMaxB = tempMaxB(c);
%                         figure
%                         clf
%                         subplot(132)
%                         hist(tempPERM,40)
%                         h = line([tempMaxB tempMaxB],[0 100]);
%                         set(h,'Color','r')
%                         
%                         tempPERM = squeeze(MinPermB(j+1,i,:));
%                         tempMinB = sort(tempPERM,'ascend');
%                         tempMinB = tempMinB(c);
%                         
%                         subplot(131)
%                         hist(tempPERM,40)
%                         h = line([tempMinB tempMinB],[0 100]);
%                         set(h,'Color','r')
%                         
%                         tempBETA = squeeze(temp(j+1,i,:));
%                         subplot(133)
%                         hist(tempBETA,40)
%                         h = line([tempMinB tempMinB],[0 1000]);
%                         set(h,'Color','r')
%                         h = line([tempMaxB tempMaxB],[0 1000]);
%                         set(h,'Color','r')
%                         
%                         tempBETA(find((tempBETA < tempMaxB) & (tempBETA > 0))) = 0;
%                         tempBETA(find((tempBETA > tempMinB) & (tempBETA < 0))) = 0;
%                         
%                         
%                         subplot(132)
%                         h = title(FileName);
%                         set(h,'Interpreter','none');
%                         
%                         I = zeros(V.dim);
%                         I(ModelInfo.Indices) = tempBETA;
%                         Vo = V;
%                         Vo.fname = fullfile(fileparts(pwd),FileName);
%                         spm_write_vol(Vo,I);
%                         
%                         
%                         
%                         
%                     end
%                     
%                 end
%             end
%         end
%     end
% end
%% PATH IMAGES

load(fullfile(fileparts(pwd),'data','AllData'));
ModelInfo.Thresholds = 0.025;
for j = 1:length(ModelInfo.Thresholds)
    c = floor(AnalysisParameters.Thresholds(j)*Nperm);
    for kk = 1:o % cycle over the number of paths
    for i = 1:m
        sMax(:,i) = sort(squeeze(MaxPermPaths(i,:,kk,:)),'descend');
        sMin(:,i) = sort(squeeze(MinPermPaths(i,:,kk,:)));
        Mx(i) = sMax(c,i);
        Mn(i) = sMin(c,i);
        % write unthresholed path estimate
        V.fname=(fullfile(fileparts(pwd),sprintf('Path%d_level%d.nii',kk,i)));
        I = zeros(V.dim);
        temp = squeeze(PointEstimate(i,1,kk,:));
         figure
        subplot(131)
        hist(sMin(:,i),40)
        h = line([Mn(i) Mn(i)],[0 100]);
        set(h,'Color','r')
                subplot(132)
        hist(sMax(:,i),40)
        h = line([Mx(i) Mx(i)],[0 100]);
        set(h,'Color','r')
        subplot(133)
        hist(temp,40)
        h = line([Mx(i) Mx(i)],[0 1000]);
        set(h,'Color','r')
          h = line([Mn(i) Mn(i)],[0 1000]);
        set(h,'Color','r')       
        
        
        I(ModelInfo.Indices) = temp;
        spm_write_vol(V,I);
        
        temp(find((PointEstimate(i,1,kk,:) < Mx(i))&(PointEstimate(i,1,kk,:) >0))) = 0;
        temp(find((PointEstimate(i,1,kk,:) > Mn(i))&(PointEstimate(i,1,kk,:) <0))) = 0;
        

        
        V.fname=(fullfile(fileparts(pwd),sprintf('Path%d_level%d_%0.4f.nii',kk,i,ModelInfo.Thresholds(j))));
        I = zeros(V.dim);
        I(ModelInfo.Indices) = temp;
        spm_write_vol(V,I);
        
    end
    end
end


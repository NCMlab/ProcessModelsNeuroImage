function WIPWriteOutImages
BaseDir = pwd;
%load(fullfile(fileparts(BaseDir),'data','AnalysisParameters'))
load(fullfile(fileparts(BaseDir),'data','AllData'))
load ../../V

AnalysisParameters = ModelInfo;
AnalysisParameters.data = [];
F = dir('Permute*');
NFiles = length(F)

load(F(1).name)
[m n o p] = size(MaxPaths);
Nperm = p*NFiles
MaxPermPaths = zeros(m,n,o,p*NFiles);
MinPermPaths = zeros(m,n,o,p*NFiles);
[mB nB pB] = size(MaxB);
MaxPermB = zeros(mB,nB,pB*NFiles);
MinPermB = zeros(mB,nB,pB*NFiles);


for i = 1:NFiles
    load(F(i).name)
    MaxPermPaths(:,:,:,(i-1)*p+1:i*p) = MaxPaths;
    MinPermPaths(:,:,:,(i-1)*p+1:i*p) = MinPaths;
    MaxPermB(:,:,(i-1)*p+1:i*p) = MaxB;
    MinPermB(:,:,(i-1)*p+1:i*p) = MinB;
end

figure(1)
hist(squeeze(MinPermPaths(1,1,1,:)),20)

figure(2)
hist(squeeze(MaxPermPaths(1,1,1,:)),20)


F = dir( 'PointEstimate*.mat')
load(F(1).name)
PointEstimate = zeros(m,n,o,length(Parameters));
for i = 1:length(Parameters)
    for kk = 1:o
        PointEstimate(:,:,kk,i) = Parameters{i}.Paths{kk};
    end
end
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
%% WRITE OUT ALL IMAGES 
Tag = 'beta';

temp = zeros([size(Parameters{1}.beta) Nvoxels]);
for i =1:Nvoxels
    temp(:,:,i) = getfield(Parameters{i},Tag);
end

for i = 1:Nvar % COLUMNS
    if sum(ModelInfo.Direct(:,i))
        %         % CONSTANT TERMS ARE ROW 1
        %         FileName = sprintf('Model%d_DEP%s_INDconst_%s.nii',i,ModelInfo.Names{i},Tag);
        %         I = zeros(V.dim);
        %         I(ModelInfo.Indices) = squeeze(temp(1,i,:));
        %         Vo = V;
        %         Vo.fname = fullfile(fileparts(pwd),FileName);
        %         spm_write_vol(Vo,I);
        for j = 1:Nvar % ROWS
            if ModelInfo.Direct(j,i)
                FileName = sprintf('Model%d_DEP%s_IND%s_%s.nii',i,ModelInfo.Names{i},ModelInfo.Names{j},Tag);
                I = zeros(V.dim);
                I(ModelInfo.Indices) = squeeze(temp(j+1,i,:));
                Vo = V;
                Vo.fname = fullfile(fileparts(pwd),FileName);
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
            I = zeros(V.dim);
            I(ModelInfo.Indices) = squeeze(temp(Nvar+2,i,:));
            Vo = V;
            Vo.fname = fullfile(fileparts(pwd),FileName);
            spm_write_vol(Vo,I);
        end
    end
end


%% WRITE OUT ALL IMAGES 
Tag = 't';

temp = zeros([size(Parameters{1}.beta) Nvoxels]);
for i =1:Nvoxels
    temp(:,:,i) = getfield(Parameters{i},Tag);
end

for i = 1:Nvar % COLUMNS
    if sum(ModelInfo.Direct(:,i))
        %         % CONSTANT TERMS ARE ROW 1
        %         FileName = sprintf('Model%d_DEP%s_INDconst_%s.nii',i,ModelInfo.Names{i},Tag);
        %         I = zeros(V.dim);
        %         I(ModelInfo.Indices) = squeeze(temp(1,i,:));
        %         Vo = V;
        %         Vo.fname = fullfile(fileparts(pwd),FileName);
        %         spm_write_vol(Vo,I);
        for j = 1:Nvar % ROWS
            if ModelInfo.Direct(j,i)
                FileName = sprintf('Model%d_DEP%s_IND%s_%s.nii',i,ModelInfo.Names{i},ModelInfo.Names{j},Tag);
                I = zeros(V.dim);
                I(ModelInfo.Indices) = squeeze(temp(j+1,i,:));
                Vo = V;
                Vo.fname = fullfile(fileparts(pwd),FileName);
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
            I = zeros(V.dim);
            I(ModelInfo.Indices) = squeeze(temp(Nvar+2,i,:));
            Vo = V;
            Vo.fname = fullfile(fileparts(pwd),FileName);
            spm_write_vol(Vo,I);
        end
    end
end











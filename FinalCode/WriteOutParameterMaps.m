function WriteOutParameterMaps(Tag,Parameters,ModelInfo,FDRFlag)
if nargin == 3
    FDRFlag = 0;
end
% create an array to hold the data of interest
% The Tag defines which parameter measure to save out.

% If a subfield is passed.
if ~isempty(findstr(Tag,'.'))
    Field = Tag(1:findstr(Tag,'.')-1);
    temp1 = getfield(Parameters{1},Field);
    temp = zeros([size(getfield(temp1,Tag(findstr(Tag,'.')+1:end))) ModelInfo.Nvoxels]);
else
    temp = zeros([size(getfield(Parameters{1},Tag)) ModelInfo.Nvoxels]);
end
% If a subfield is passed.
for i = 1:ModelInfo.Nvoxels
    if ~isempty(findstr(Tag,'.'))
        Field = Tag(1:findstr(Tag,'.')-1);
        temp1 = getfield(Parameters{i},Field);
        temp(:,:,i) = getfield(temp1,Tag(findstr(Tag,'.')+1:end));
    else
        temp(:,:,i) = getfield(Parameters{i},Tag);
    end
end

if strcmp(Tag,'BCaCI.p') && FDRFlag
    Tag = sprintf('%s_%s',Tag,'FDR');
    % Only use the first threshold
    kk = 1;
    alpha = ModelInfo.Thresholds(kk);
    % Cycle over each parameter
    for i = 1:ModelInfo.Nvar
        if sum(ModelInfo.Direct(:,i))
            % CONSTANT TERMS ARE ROW 1
            % cycle over the rows in the model
            for j = 1:ModelInfo.Nvar
                if ModelInfo.Direct(j,i)
                    I = sort(squeeze(temp(j+1,i,:)));
                    [FDRp, FDRpN] = FDR(I,alpha);
                    if isempty(FDRp)
                        FDRp = 1/ModelInfo.Nvoxels*alpha;
                    end
                    temp(j+1,i,find(temp(j+1,i,:) > FDRp)) = 0;
                end
            end
        end
    end
end

% cycle over the columns in the model
for i = 1:ModelInfo.Nvar 
    if sum(ModelInfo.Direct(:,i))
        % CONSTANT TERMS ARE ROW 1
        % cycle over the rows in the model
        for j = 1:ModelInfo.Nvar 
            if ModelInfo.Direct(j,i)
                % create the filename describing the dependent and
                % independent effects
                FileName = sprintf('Model%d_DEP%s_IND%s_%s.nii',i,ModelInfo.Names{i},ModelInfo.Names{j},Tag);
                
                I = zeros(ModelInfo.DataHeader.dim);
                I(ModelInfo.Indices) = squeeze(temp(j+1,i,:));
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
            I = zeros(ModelInfo.DataHeader.dim);
            I(ModelInfo.Indices) = squeeze(temp(ModelInfo.Nvar+2,i,:));
            % Create the header for this image
            Vo = ModelInfo.DataHeader;
            Vo.fname = fullfile(ModelInfo.ResultsPath,FileName);
            spm_write_vol(Vo,I);
        end
    end
end

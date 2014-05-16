function WriteOutParameterMaps(Tag,Parameters,ModelInfo)

% create an array to hold the data of interest
% The Tag defines which parameter measure to save out.

temp = zeros([size(getfield(Parameters{1},Tag)) ModelInfo.Nvoxels]);
for i =1:ModelInfo.Nvoxels
    temp(:,:,i) = getfield(Parameters{i},Tag);
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

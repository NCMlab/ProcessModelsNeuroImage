function WriteOutParameterMaps(Tag)


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

function WriteOutSurfaceParameterMaps(Tag,Parameters,ModelInfo,FDRFlag)


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
    % cycle over each entry in the data and apply FDR to each
    [m n p] = size(temp);
    % Cycle over each parameter
    for i = 1:n%ModelInfo.Nvar
        %if sum(ModelInfo.Direct(:,i))
        % CONSTANT TERMS ARE ROW 1
        % cycle over the rows in the model
        for j = 1:m%ModelInfo.Nvar
            %if ModelInfo.Direct(j,i)
            
            I = sort(squeeze(temp(j,i,:)));
            if sum(abs(I)) > 0
                [FDRp, FDRpN] = FDR(I,alpha);
                if isempty(FDRp)
                    FDRp = 1/ModelInfo.Nvoxels*alpha;
                end
                temp(j,i,find(temp(j,i,:) > FDRp)) = 0;
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
                FileName = sprintf('Model%d_DEP%s_IND%s_%s',i,ModelInfo.Names{i},ModelInfo.Names{j},Tag);
                % Working with surface maps
                Ilh = zeros(length(ModelInfo.lhVertices),1);
                Irh = zeros(length(ModelInfo.rhVertices),1);
                LHindices = 1:length(ModelInfo.lhVertices);
                RHindices = 1:length(ModelInfo.rhVertices);
                Ilh(LHindices) = squeeze(temp(j+1,i,LHindices));
                Irh(RHindices) = squeeze(temp(j+1,i,length(LHindices)+RHindices));
                % Write data to files
                WriteVertexFiles(Ilh,fullfile(ModelInfo.ResultsPath,[FileName '.lh.asc']),ModelInfo.lhVertices')
                WriteVertexFiles(Irh,fullfile(ModelInfo.ResultsPath,[FileName '.rh.asc']),ModelInfo.rhVertices')
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
                % Working with surface maps
                Ilh = zeros(length(ModelInfo.lhVertices),1);
                Irh = zeros(length(ModelInfo.rhVertices),1);
                LHindices = 1:length(ModelInfo.lhVertices);
                RHindices = 1:length(ModelInfo.rhVertices);
                Ilh(LHindices) = squeeze(temp(j+1,i,LHindices));
                Irh(RHindices) = squeeze(temp(j+1,i,length(LHindices)+RHindices));
                % Write data to files
                WriteVertexFiles(Ilh,fullfile(ModelInfo.ResultsPath,[FileName '.lh.asc']),ModelInfo.lhVertices)
                WriteVertexFiles(Irh,fullfile(ModelInfo.ResultsPath,[FileName '.rh.asc']),ModelInfo.rhVertices)

        
        end
    end
end


function WriteVertexFiles(Data,File,Vertices)
fid = fopen(File,'w');
if fid < 3
    error('Cannot open file!');
else
    for i = 1:length(Data)
        fprintf(fid,'%03d %0.5f %0.5f %0.5f %0.5f\n',Vertices(i,1),Vertices(i,2),Vertices(i,3),Vertices(i,4),Data(i,1));
    end
end
fclose(fid);
fprintf(1,'Data written to: %s\n',File);



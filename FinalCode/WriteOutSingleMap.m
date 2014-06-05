function WriteOutSingleMap(Tag,Parameters,ModelInfo,FDRFlag)
if nargin == 3
    FDRFlag = 0;
end



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

if strcmp(Tag,'BCaCI.PathsP') && FDRFlag
    Tag = sprintf('%s_%s',Tag,'FDR');
    % Only use the first threshold
    kk = 1;
    alpha = ModelInfo.Thresholds(kk);
    % Cycle over each parameter
    I = sort(squeeze(temp));
    [FDRp, FDRpN] = FDR(I,alpha);
    if isempty(FDRp)
        FDRp = 1/ModelInfo.Nvoxels*alpha;
    end
    temp(find(temp) > FDRp) = 0;
end

FileName = fullfile(ModelInfo.ResultsPath,sprintf('Path_%s.nii',Tag));

I = zeros(ModelInfo.DataHeader.dim);
I(ModelInfo.Indices) = squeeze(temp);
% Create the header for this image
Vo = ModelInfo.DataHeader;
Vo.fname = FileName;
spm_write_vol(Vo,I);


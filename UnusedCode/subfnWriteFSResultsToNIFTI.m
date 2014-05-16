function subfnWriteFSResultsToNIFTI(OutData,OutFolder,MeasureOfInterest,Header)

CodeFolder = '/Users/jason/Dropbox/SteffenerColumbia/Scripts/ProcessModelsNeuroImage/FreeSurferFiles';

CodeFolder = '/share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FreeSurferFiles';
% Use the Freesurfer lookup table to map location names to spots in the
% brain. The file has been cleaned so that it does not have any of the
% descriptive breaks that are in the original
D = textread(fullfile(CodeFolder,'FreeSurferColorLUTclean.txt'),'%s','delimiter',' ');

NCol = 6;
NRow = length(D)/NCol;
Data = reshape(D,NCol,NRow)';
label = str2num(char(Data(:,1)));
% create an array of all the names of the labeled regions
name = Data(:,2);
% From the names strip off the CTX labels
for i = 1:length(name)
    if strfind(name{i},'ctx') == 1 
        name{i} = strrep(name{i},'ctx-','');
    end
end


% convert dashes to underscores
for i = 1:length(name)
    temp = name{i};
    fDash = findstr(temp,'-');
    if ~isempty(fDash)
        temp(fDash) = '_';
        name{i} = temp;
    end
end

for i = 1:length(Header)
    temp = Header{i};
    fDash = findstr(temp,'-');
    if ~isempty(fDash)
        temp(fDash) = '_';
        Header{i} = temp;
    end
end

    
        
% Load up the example images
% These are taken from a single subject and they really shuld be rotated
% correctly. 
%P = fullfile('/share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FreeSurferFiles','raparc+asegROT.nii');
P = fullfile(CodeFolder,'STANDARD_aparc+aseg.nii');
V = spm_vol(P);
I = spm_read_vols(V);
% for all locations in the Freesurfer header file, which comes from
% aparcstats2table, map these to values in the image (label)
% for i = 1:length(Header)
%     fprintf(1,'%d\t%s\n',i,Header{i});
% end

for i = 1:length(name)
    fprintf(1,'%d\t%s\n',i,name{i});
end

FSNameAndLabel = {};
for i = 1:length(Header)
    for j = 1:length(name)
        if strfind(Header{i},[name{j}])
            FSNameAndLabel{i}.name = [Header{i}];
            FSNameAndLabel{i}.label = label(j);
            F = find(I == FSNameAndLabel{i}.label);
            FSNameAndLabel{i}.voxels = F;
            break
        else
            FSNameAndLabel{i} = [];
        end
        
    end
end


List = {};
index = 1;
for i = 1:length(FSNameAndLabel)
    if ~isempty(FSNameAndLabel{i})
        List{index} = FSNameAndLabel{i}.name;
        index = index + 1;
    end
end






% Create a mapping of the brain regions back to the parameters in the
% results.
%OutData = {};
index = length(OutData) + 1;
OutData{index}.name = 'BrainToParameters';
OutData{index}.data = ones(84,1);
OutData{index}.data = data;
%[1:length(OutData{index-1}.data)];
OutData{index}.dataType = 4;

% Write the images out
% Check for the thresholdindex values and the order values. These will be
% used to create 4-D images of the probe values.
FilesToRemove = {};
Str = '';
for i = 1:length(OutData)
%     thresholdIndex = 0;
    order = 0;
    % check for threshold index
%     if isfield(OutData{i},'thresholdIndex')
%         thresholdIndex = OutData{i}.thresholdIndex;
%     end
    if isfield(OutData{i},'order')
        order = OutData{i}.order;
    end
    
    if order == 1
        Str
        unix(Str);
        Str = sprintf('fslmerge -t %s_%s_4D',MeasureOfInterest,OutData{i}.name);
    end
    
    Vo = V;
    Vo.dt = [16 0];    
    Vo.fname = fullfile(OutFolder,[OutData{i}.name '_' MeasureOfInterest '.nii']);
    Y = zeros(Vo.dim);
    for j = 1:length(FSNameAndLabel)
        if ~isempty(FSNameAndLabel{j})
            Y(FSNameAndLabel{j}.voxels) = OutData{i}.data(j);
        end
    end
    spm_write_vol(Vo,Y);
    if order > 0
        Str = sprintf('%s %s',Str,Vo.fname);
        FilesToRemove{length(FilesToRemove)+1} = Vo.fname;
    end
end

unix(Str);
for i = 1:length(FilesToRemove)
    Str = sprintf('rm %s',FilesToRemove{i});
    unix(Str);
end
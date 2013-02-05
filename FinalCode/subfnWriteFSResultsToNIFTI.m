function subfnWriteFSResultsToNIFTI(OutData,OutFolder,BasePath)

FSPath = fullfile(BasePath,'FreeSurferFiles');
% Load up the Header file used with the data
load(fullfile(BasePath,'FSheader'))
% Use the Freesurfer lookup table to map location names to spots in the
% brain. The file has been cleaned so that it does not have any of the
% descriptive breaks that are in the original
D = textread(fullfile(BasePath,'FreeSurferColorLUTclean.txt'),'%s','delimiter',' ');

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
        
    
        
% Load up the example images
% These are taken from a single subject and they really shuld be rotated
% correctly. 
P = fullfile(BasePath,'raparc+asegROT.nii');
V = spm_vol(P);
I = spm_read_vols(V);
% for all locations in the Freesurfer header file, which comes from
% aparcstats2table, map these to values in the image (label)
% for i = 1:length(Header)
%     fprintf(1,'%d\t%s\n',i,Header{i});
% end

MeasureOfInterest = 'thickness';


FSNameAndLabel = {};
for i = 1:length(Header)
    for j = 1:length(name)
        if strmatch(Header{i},[name{j} '_' MeasureOfInterest])
            FSNameAndLabel{i}.name = [name{j} '_' MeasureOfInterest];
            FSNameAndLabel{i}.label = label(j);
            F = find(I == FSNameAndLabel{i}.label);
            FSNameAndLabel{i}.voxels = F;
            break
        else
             FSNameAndLabel{i} = [];
        end
        
    end
end

% Write the images out
% To Do: make sure to change the output name and folders
for i = 1:length(OutData)
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
end

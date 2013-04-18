function [outNames Hemi MeasurementType] = subfnModifyFreesurferNames(inNames)

N = length(inNames);
Hemi = cell(N,1);
outNames = cell(N,1);
MeasurementType = cell(N,1);
FS = createNames;

for i = 1:N
    if ~isempty(findstr(inNames{i},'lh'))
        Hemi{i} = 'L';
        % strip off the hemispher part in the front 
        % strip off the thickness part at the end
        outNames{i} = FindCorrectName(FS,inNames{i}(4:end-length('thickness')-1));
        MeasurementType{i} = 'T';
    elseif ~isempty(findstr(inNames{i},'rh'))
        Hemi{i} = 'R';
        outNames{i} = FindCorrectName(FS,inNames{i}(4:end-length('thickness')-1));
        MeasurementType{i} = 'T';
    elseif ~isempty(findstr(inNames{i},'Left'))
        Hemi{i} = 'L';
        outNames{i} = FindCorrectName(FS,inNames{i}(6:end));
        MeasurementType{i} = 'V';
    elseif ~isempty(findstr(inNames{i},'Right'))
        Hemi{i} = 'R';
        outNames{i} = FindCorrectName(FS,inNames{i}(7:end));
        MeasurementType{i} = 'V';
    end
end


function outName = FindCorrectName(FS,inName)
for i = 1:length(FS)
    if strmatch(inName,FS{i,1}.in);
        if isfield(FS{i,1},'out')
            outName = FS{i,1}.out;
        else
            outName = FS{i,1}.in;
        end
        break
    end
end

function FS = createNames
FS{1,1}.in = 'Accumbens-area';
FS{1,1}.out = 'Accumbens';
FS{2,1}.in =     'Amygdala';
FS{3,1}.in =     'Caudate';
FS{4,1}.in =     'Hippocampus';
FS{5,1}.in =     'Pallidum';
FS{6,1}.in =     'Putamen';
FS{7,1}.in =     'Thalamus-Proper';
FS{7,1}.out =     'Thalamus';
FS{8,1}.in =     'VentralDC';
FS{9,1}.in =     'bankssts';
FS{10,1}.in =     'caudalanteriorcingulate';
FS{10,1}.out =     'Ant. Caudal Cingulate';
FS{11,1}.in =     'caudalmiddlefrontal';
FS{11,1}.out =     'Mid. Caudal Frontal';
FS{12,1}.in =     'cuneus';
FS{12,1}.out =     'Cuneus';
FS{13,1}.in =     'entorhinal';
FS{13,1}.out =     'Entorhinal';
FS{14,1}.in =    'frontalpole';
FS{14,1}.out =    'Frontal Pole';
FS{15,1}.in =    'fusiform';
FS{15,1}.out =    'Fusiform';
FS{16,1}.in =     'inferiorparietal';
FS{16,1}.out =     'Inf. Parietal';
FS{17,1}.in =     'inferiortemporal';
FS{17,1}.out =     'Inf. Temporal';
FS{18,1}.in =     'insula';
FS{18,1}.out =     'Insula';
FS{19,1}.in =     'isthmuscingulate';
FS{19,1}.out =     'Isthmus of Cingulate';
FS{20,1}.in =     'lateraloccipital';
FS{20,1}.out =     'Lat. Occipital';
FS{21,1}.in =     'lateralorbitofrontal';
FS{21,1}.out =     'Lat. Orbitofrontal';
FS{22,1}.in =     'lingual';
FS{22,1}.out =     'Lingual';
FS{23,1}.in =     'medialorbitofrontal';
FS{23,1}.out =     'Med. Orbitofrontal';
FS{24,1}.in =     'middletemporal';
FS{24,1}.out =     'Mid. Temporal';
FS{25,1}.in =     'paracentral';
FS{25,1}.out =     'Paracentral';
FS{26,1}.in =     'parahippocampal';
FS{26,1}.out =     'Parahippocampal';
FS{27,1}.in =     'parsopercularis';
FS{27,1}.out =     'Parsopercularis';
FS{28,1}.in =     'parsorbitalis';
FS{28,1}.out =     'Parsorbitalis';
FS{29,1}.in =     'parstriangularis';
FS{29,1}.out =     'Parstriangularis';
FS{30,1}.in =     'pericalcarine';
FS{30,1}.out =     'Pericalcarine';
FS{31,1}.in =     'postcentral';
FS{31,1}.out =     'Post. Central';
FS{32,1}.in =     'posteriorcingulate';
FS{32,1}.out =     'Post. Cingulate';
FS{33,1}.in =     'precentral';
FS{33,1}.out =     'Precentral';
FS{34,1}.in =     'precuneus';
FS{34,1}.out =     'Precuneus';
FS{35,1}.in =     'rostralanteriorcingulate';
FS{35,1}.out =     'Rost. Ant. Cingulate';
FS{36,1}.in =     'rostralmiddlefrontal';
FS{36,1}.out =     'Rost. Mid. Frontal';
FS{37,1}.in =     'superiorfrontal';
FS{37,1}.out =     'Sup. Frontal';
FS{38,1}.in =     'superiorparietal';
FS{38,1}.out =     'Sup. Parietal';
FS{39,1}.in =     'superiortemporal';
FS{39,1}.out =     'Sup. Temporal';
FS{40,1}.in =     'supramarginal';
FS{40,1}.out =     'Supramarginal';
FS{41,1}.in =     'temporalpole';
FS{41,1}.out =     'Temp. Pole';
FS{42,1}.in =     'transversetemporal';
FS{42,1}.out =     'Trans. Temporal';
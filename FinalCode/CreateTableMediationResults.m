function CreateTableMediationResults(LocalMaxima,InputImage)
% Input is a text file created from FSL
% InputImage is the file used to create the list of local maxima


[PathName FileName] = fileparts(InputImage);
if ~isempty(PathName)
    cd(PathName);
end
load(fullfile(PathName,'data','ModelInfo'));

Vpath = spm_vol(InputImage);
Ipath = spm_read_vols(Vpath);
% Find which path this is and the point estimate of it
FUnders = findstr(FileName,'_');
WhichPath = FileName(1:FUnders(1)-1);
% Find the Parameter Estimate
PathPE = dir([WhichPath '_PointEstimate.nii']);
VpathPE = spm_vol(PathPE.name);
IPathPE = spm_read_vols(VpathPE);
% Read in the local maxima file and all locations
[ID, Value, x, y, z] = textread(LocalMaxima,'%d %f %d %d %d','headerlines',1);
Nloc = length(ID);
%% Find all the parameter images for this result
% cycle over a large number of models
AllValuesBeta = [];
AllValuesT = [];
AllValuesPath = [];
DEPNames = {};
INDNames = {};
fid = 1;
for i = 1:10
    F = dir(fullfile(PathName,sprintf('Model%d_*_beta.nii',i)));
    for j = 1:length(F)
        FindUnder = findstr(F(j).name,'_');
        FindDEP = findstr(F(j).name,'DEP');
        FindIND = findstr(F(j).name,'IND');
        % Create list of column headings
        DEPNames{length(DEPNames)+1} = F(j).name(FindDEP+3:FindUnder(min(find(FindUnder>FindDEP)))-1);
        INDNames{length(INDNames)+1} = F(j).name(FindIND+3:FindUnder(min(find(FindUnder>FindIND)))-1);
        Ibeta = spm_read_vols(spm_vol(fullfile(PathName,F(j).name)));
        %BCaCI.Z
        Pt = fullfile(PathName,[strrep(F(j).name,'_beta.','_t.')]);
        %Pt = fullfile(PathName,[strrep(F(j).name,'beta','BCaCI.Z')]);
        Vt = spm_vol(Pt);
        It = spm_read_vols(Vt);
        values = zeros(Nloc,1);
        for k = 1:Nloc
            XYZ = inv(Vt.mat)*[[x(k) y(k) z(k) 1]]';
            valuesB(k) = Ibeta(XYZ(1),XYZ(2),XYZ(3));
            valuesT(k) = It(XYZ(1),XYZ(2),XYZ(3));
        end
        AllValuesBeta(:,size(AllValuesBeta,2)+1) = valuesB;
        AllValuesT(:,size(AllValuesT,2)+1) = valuesT;
    end
end
AllValuesPathP = [];
AllValuesPathPE = [];

for k = 1:Nloc
    XYZ = inv(Vt.mat)*[[x(k) y(k) z(k) 1]]';
    AllValuesPathP(k) = Ipath(XYZ(1),XYZ(2),XYZ(3));
    AllValuesPathPE(k) = IPathPE(XYZ(1),XYZ(2),XYZ(3));
end
%% Write dependent names
[PathName FileName] = fileparts(LocalMaxima);
OutFile = fullfile(pwd,['MedDetails_' FileName '.csv']);
fid = fopen(OutFile,'w');
if fid < 0
    errordlg('Cannot create new file, printing to screen.');
    fid = 1;
end
% write top row
fprintf(fid,'Input Image, %s\n',fullfile(pwd,InputImage));
fprintf(fid,',,,,,,');
for j = 1:size(DEPNames,2)
    fprintf(fid,'beta,');
end
for j = 1:size(DEPNames,2)
    fprintf(fid,'t,');
end
fprintf(fid,'\n');
% write second row
fprintf(fid,',,,,,,');
for j = 1:size(DEPNames,2)
    fprintf(fid,'%s,',DEPNames{j});
end
for j = 1:size(DEPNames,2)
    fprintf(fid,'%s,',DEPNames{j});
end
fprintf(fid,'\n');
% write third row
fprintf(fid,'ClusterIndex,PathPE,Path_Pvalue,x,y,z,');
for j = 1:size(INDNames,2)
    fprintf(fid,'%s,',INDNames{j});
end
for j = 1:size(INDNames,2)
    fprintf(fid,'%s,',INDNames{j});
end
fprintf(fid,'\n');
% Write out data

for i = 1:Nloc
    fprintf(fid,'%d,',ID(i));
    fprintf(fid,'%0.6f,',AllValuesPathPE(i));
    fprintf(fid,'%0.4f,',AllValuesPathP(i));
    fprintf(fid,'%d,',x(i));
    fprintf(fid,'%d,',y(i));
    fprintf(fid,'%d,',z(i));
    for j = 1:size(INDNames,2)
        fprintf(fid,'%0.6f,',AllValuesBeta(i,j));
    end
    for j = 1:size(INDNames,2)
        fprintf(fid,'%0.6f,',AllValuesT(i,j));
    end

    fprintf(fid,'\n');
end
if fid > 1
    fclose(fid);
    fprintf(1,'Results written to file: \n \t%s\n',OutFile);
end
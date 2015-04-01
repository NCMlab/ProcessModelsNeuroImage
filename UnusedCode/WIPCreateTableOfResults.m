% Select PATH image
P = spm_select(1,'image');
[PathName FileName] = fileparts(P);
cd(PathName);



Vpath = spm_vol(P);
%% Negative
HeightThreshold = 0;
ExtentThreshold = 10;
[xSPM hReg] = SPMDisplayThresholdImage(P,HeightThreshold,ExtentThreshold);
[XYZmm] = FindAALandBAForSPMResults(xSPM, hReg,[FileName]);
Nloc = size(XYZmm,1);


% Find all the parameter images for this result
% cycle over a large number of models
AllValuesBeta = [];
AllValuesT = [];
DEPNames = {};
INDNames = {};
for i = 1:10
    F = dir(fullfile(PathName,sprintf('Model%d_*_B.nii',i)));
    for j = 1:length(F)
        FindUnder = findstr(F(j).name,'_');
        FindDEP = findstr(F(j).name,'DEP');
        FindIND = findstr(F(j).name,'IND');
        DEPNames{length(DEPNames)+1} = F(j).name(FindDEP+3:FindUnder(min(find(FindUnder>FindDEP)))-1);
        INDNames{length(INDNames)+1} = F(j).name(FindIND+3:FindUnder(min(find(FindUnder>FindIND)))-1);
        Ibeta = spm_read_vols(spm_vol(fullfile(PathName,F(j).name)));
        Pt = fullfile(PathName,[strrep(F(j).name,'beta','t')]);
        Vt = spm_vol(Pt);
        It = spm_read_vols(Vt);
        values = zeros(Nloc,1);
        for k = 1:Nloc
            XYZ = inv(Vt.mat)*[XYZmm(k,:) 1]';
            valuesB(k) = Ibeta(XYZ(1),XYZ(2),XYZ(3));
            valuesT(k) = It(XYZ(1),XYZ(2),XYZ(3));
        end
        AllValuesBeta(:,size(AllValuesBeta,2)+1) = valuesB;
        AllValuesT(:,size(AllValuesT,2)+1) = valuesT;
    end
end
Tthresh = 1.96;
% Write dependent names
fid = 1;
for j = 1:size(DEPNames,2)
    fprintf(fid,'%s\t',DEPNames{j});
end
fprintf(fid,'\n');
for j = 1:size(INDNames,2)
    fprintf(fid,'%s\t',INDNames{j});
end
fprintf(fid,'\n');
for i = 1:Nloc
    for j = 1:size(INDNames,2)
        
        fprintf(fid,'%0.2f',AllValuesBeta(i,j));
        if abs(AllValuesT(i,j))>Tthresh
            fprintf(fid,'*\t');
        else
            fprintf(fid,'\t');
        end
    end
    fprintf(fid,'\n');
end

%% Positive
HeightThreshold = 10e-6;

[xSPM hReg] = SPMDisplayThresholdImage(P,HeightThreshold,ExtentThreshold);
[XYZmm] = FindAALandBAForSPMResults(xSPM, hReg,[FileName]);
Nloc = size(XYZmm,1);


% Find all the parameter images for this result
% cycle over a large number of models
AllValuesBeta = [];
AllValuesT = [];
DEPNames = {};
INDNames = {};
for i = 1:10
    F = dir(fullfile(PathName,sprintf('Model%d_*_B.nii',i)));
    for j = 1:length(F)
        FindUnder = findstr(F(j).name,'_');
        FindDEP = findstr(F(j).name,'DEP');
        FindIND = findstr(F(j).name,'IND');
        DEPNames{length(DEPNames)+1} = F(j).name(FindDEP+3:FindUnder(min(find(FindUnder>FindDEP)))-1);
        INDNames{length(INDNames)+1} = F(j).name(FindIND+3:FindUnder(min(find(FindUnder>FindIND)))-1);
        Ibeta = spm_read_vols(spm_vol(fullfile(PathName,F(j).name)));
        Pt = fullfile(PathName,[strrep(F(j).name,'beta','t')]);
        Vt = spm_vol(Pt);
        It = spm_read_vols(Vt);
        values = zeros(Nloc,1);
        for k = 1:Nloc
            XYZ = inv(Vt.mat)*[XYZmm(k,:) 1]';
            valuesB(k) = Ibeta(XYZ(1),XYZ(2),XYZ(3));
            valuesT(k) = It(XYZ(1),XYZ(2),XYZ(3));
        end
        AllValuesBeta(:,size(AllValuesBeta,2)+1) = valuesB;
        AllValuesT(:,size(AllValuesT,2)+1) = valuesT;
    end
end
Tthresh = 1.96;
% Write dependent names
fid = 1;
for j = 1:size(DEPNames,2)
    fprintf(fid,'%s\t',DEPNames{j});
end
fprintf(fid,'\n');
for j = 1:size(INDNames,2)
    fprintf(fid,'%s\t',INDNames{j});
end
fprintf(fid,'\n');
for i = 1:Nloc
    for j = 1:size(INDNames,2)
        
        fprintf(fid,'%0.2f',AllValuesBeta(i,j));
        if abs(AllValuesT(i,j))>Tthresh
            fprintf(fid,'*\t');
        else
            fprintf(fid,'\t');
        end
    end
    fprintf(fid,'\n');
end


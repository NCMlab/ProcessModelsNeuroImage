function OutFileName = subfnWriteDatatoCSV(Location)


BasePath = '/share/users/js2746_Jason/Studies/ApoEStudy'

load(fullfile(BasePath,'FSheader.mat'));
for i = 1:length(Header)
    fprintf(1,'%d\t%s\n',i,Header{i});
end

SelectedPath = spm_select(1,'dir');
cd(SelectedPath)
if exist('data_0001.mat')
    load data_0001
else
    errordlg('this folder does not have the required AnalysticParameters.mat file');
end


OutFileName = fullfile(SelectedPath,['Location_' Header{Location} '.csv']);
fid = fopen(OutFileName,'w');

NSub = length(data.X);
if ~isempty(data.names.X)
    fprintf(fid,'%s,',data.names.X);
end
if ~isempty(data.names.M)
    for i = 1:length(data.names.M)
        fprintf(fid,'%s,',data.names.M{i});
    end
end
if ~isempty(data.names.Y)
    fprintf(fid,'%s,',data.names.Y);
end
if ~isempty(data.names.V)
    fprintf(fid,'%s,',data.names.V);
end
if ~isempty(data.names.W)
    fprintf(fid,'%s,',data.names.W);
end
if  ~isempty(data.names.COV)
    for i = 1:length(data.names.COV)
        fprintf(fid,'%s,',data.names.COV{i});
    end
end
if ~isempty(data.names.Q)
    fprintf(fid,'%s,',data.names.Q);
end
if ~isempty(data.names.R)
    fprintf(fid,'%s,',data.names.R);
end
fprintf(fid,'\n');
for i = 1:NSub
    if ~isempty(data.names.X)
        fprintf(fid,'%0.4f,',data.X(i));
    end
    if ~isempty(data.names.M)
        for j = 1:length(data.names.M)
            fprintf(fid,'%0.4f,',data.M(i,j,Location));
        end
    end
    if ~isempty(data.names.Y)
        fprintf(fid,'%0.4f,',data.Y(i));
    end
    if ~isempty(data.names.V)
        fprintf(fid,'%0.4f,',data.V(i));
    end
    if ~isempty(data.names.W)
        fprintf(fid,'%0.4f,',data.W(i));
    end
    if  ~isempty(data.names.COV)
        for j = 1:length(data.names.COV)
            fprintf(fid,'%0.4f,',data.COV(i,j));
        end
    end
    if ~isempty(data.names.Q)
        fprintf(fid,'%0.4f,',data.Q(i));
    end
    if ~isempty(data.names.R)
        fprintf(fid,'%0.4f,',data.R(i));
    end
    fprintf(fid,'\n');
end
fclose(fid);

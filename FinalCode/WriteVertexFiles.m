function WriteVertexFiles(Data,File,Vertices)
fid = fopen(File,'w');
if fid < 3
    error('Cannot open file!');
else
    for i = 1:length(Data)
        fprintf(fid,'%03d %0.5f %0.5f %0.5f %0.5f\n',Vertices(1,i),Vertices(2,i),Vertices(3,i),Vertices(4,i),Data(i,1));
    end
end
fclose(fid)
end
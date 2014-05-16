function Zdata = subfnStandardizeData(data)
Zdata = data;
Zdata.Y = standardize(data.Y);
for i = 1:size(data.COV,2)
    Zdata.COV(:,i) = standardize(data.COV(:,i));
end
for i = 1:size(data.M,2)
    for j = 1:size(data.M,3)
        Zdata.M(:,i,j) = standardize(data.M(:,i,j));
    end
end
ZSdata.Y = standardize(data.Y);
if numel(unique(data.X)) == 2
    Zdata.X = data.X;
else
    Zdata.X = standardize(data.X);
end

function Zx = standardize(x)
Zx = (x - mean(x))/std(x);
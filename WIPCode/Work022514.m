
PointEstResults = WIPsubfnVoxelWiseProcessBatch(AllData)


Data1 = AllData;
Data1.M = Data1.M(:,:,1:1000);
Data1.Indices = Data1.Indices(1:1000);
PointEstResults = WIPsubfnVoxelWiseProcessBatch(Data1)
Nvox = size(Data1.M,3);




Parameters2 = WIPsubfnVoxelWiseProcessBatch(Data1,1);

Nboot = 20;
Parameters
FieldNames = {'beta' 'B' 'Paths'};

PointEstResults = WIPsubfnVoxelWiseProcessBatch(Data1)

BootStrap = {};
for i = 1:length(FieldNames)
    Value = getfield(PointEstResults{1},FieldNames{i});
    if iscell(Value)
        BlankValue = cell(size(Value),Nvox,Nboot);
    else
        BlankValue = zeros([size(Value) Nvox Nboot]);
    end
    BootStrap = setfield(BootStrap,FieldNames{i},BlankValue);
end


for i = 1:Nboot
    BootParameters = WIPsubfnVoxelWiseProcessBatch(Data1,1);
    for j = 1:Nvox
        BootStrap.beta(:,:,j,i) = BootParameters{j}.beta;
        BootStrap.B(:,:,j,i) = BootParameters{j}.B;
        BootStrap.Paths{:,j,i} = BootParameters{j}.Paths{1};
    end
end

PE = zeros(Nvox,1);
BB = zeros(Nvox,Nboot);
for j = 1:Nvox
    PE(j) = PointEstResults{j}.Paths{1}(1);
    
    for k = 1:Nboot
        BB(j,k) = BootStrap.Paths{1,j,k}(1);
    end
end


maxBB = max(BB)';
SmaxBB = sort(maxBB);
a = 0.05;
SmaxBB(floor(a*Nboot))
F = find(PE>SmaxBB(floor(a*Nboot)))


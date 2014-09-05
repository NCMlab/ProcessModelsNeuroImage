function BoLBBCaCI = BagOfBootStrapFunction(data,Nboot,FieldNames,JackKnife,lambda)
fprintf(1,'Calculating the Bag of Little Bootstrap.\n')
PointEstResults  = FitProcessModel(data);

% initialize bootstrap values based on the size of the results structure.
% FieldNames to bootstrap
Nsubsamples = 500;

%%
BootStrap = {};
for i = 1:length(FieldNames)
    Value = getfield(PointEstResults,FieldNames{i});
    if iscell(Value)
        BlankValue = cell(size(Value,1),Nboot);
    else
        BlankValue = zeros([size(Value) Nboot]);
    end
    BootStrap = setfield(BootStrap,FieldNames{i},BlankValue);
end
%%

N = data.Nsub;
Nsubsample = ceil(N*lambda);
BagBootStrapCI = cell(Nsubsamples,1);
count = 1;

while count < Nsubsamples+1
    try
        SubsampleIndex = randperm(N);
        SubsampleIndex = SubsampleIndex(1:Nsubsample);
        SubsampleData = data;
        SubsampleData.data = SubsampleData.data(SubsampleIndex,:);
        SubSampleBootStrap = BootStrap;
        for j = 1:Nboot
            temp = SubsampleData;
            % create the bootstrap resmples having the length of the original
            % data, but only using the sample indices drawn from the subsampled
            % data set.
            BootSamp =  floor(Nsubsample*rand(data.Nsub,1))+1;
            temp.data = temp.data(BootSamp,:);
            Results = FitProcessModel(temp);
            SubSampleBootStrap.beta(:,:,j) = Results.beta;
            SubSampleBootStrap.B(:,:,j) = Results.B;
            SubSampleBootStrap.Paths{j} = Results.Paths;
        end
        % Now calculate the confidence intervals for this subsample
        BCaCI = CreateBCaCI(PointEstResults,SubSampleBootStrap,JackKnife,data.Thresholds);
        BagBootStrapCI{count} = BCaCI;
        count = count + 1;
        %fprintf(1,'%d\n',count);
    catch me
        
    end
end
% Average across the confidence intervals

Sumbeta = zeros(size(BagBootStrapCI{1}.beta));
SumB = zeros(size(BagBootStrapCI{1}.B));
SumPath = zeros(size(BagBootStrapCI{1}.Paths));
for i = 1:Nsubsamples
    Sumbeta = Sumbeta + BagBootStrapCI{i}.beta;
    SumB = SumB + BagBootStrapCI{i}.B;
    SumPath = SumPath + BagBootStrapCI{i}.Paths;
end
BoLBBCaCI.beta = Sumbeta./Nsubsamples;
BoLBBCaCI.B = SumB./Nsubsamples;
BoLBBCaCI.Path = SumPath./Nsubsamples;


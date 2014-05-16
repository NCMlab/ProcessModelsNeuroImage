function BootStrap = BootStrapFunction(data,Nboot,FieldNames)

PointEstResults  = FitProcessModel(data);

% initialize bootstrap values based on the size of the results structure.
% FieldNames to bootstrap

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
% find the number of groups if stratified resampling will be done
if ~isempty(data.STRAT)
    Gr1 = find(data.STRAT == 0);
    NGr1 = length(Gr1);
    Gr2 = find(data.STRAT == 1);
    NGr2 = length(Gr2);
else
    NGr1 = [];
    NGr2 = [];
end
% start the bootstrap loop using parallel processing
for i = 1:Nboot
    % this is needed for the parallel processing to work
    temp = data;
    % create the resamples
    if isempty(temp.STRAT)
        Samp =  floor(data.Nsub*rand(data.Nsub,1))+1;
    else
        Samp1 = floor(NGr1*rand(NGr1,1))+1;
        Samp2 = floor(NGr2*rand(NGr2,1))+1+NGr1;
        Samp = [Samp1; Samp2];
    end
    temp.data = temp.data(Samp,:);
    Results = FitProcessModel(temp);
    BootStrap.beta(:,:,i) = Results.beta;
    BootStrap.B(:,:,i) = Results.B;
    BootStrap.Paths{i} = Results.Paths;
end

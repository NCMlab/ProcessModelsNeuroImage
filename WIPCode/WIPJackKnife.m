function JackKnife = WIPJackKnife(Model,Results)
tempModel = Model;
tempModel.N = Model.N - 1;
FieldNames = {'beta' 'B' 'Paths'};

% create a structure to hold the jack-knife calculations
JackKnife = {};
for i = 1:length(FieldNames)
    Value = getfield(Results,FieldNames{i});
    if iscell(Value)
        BlankValue = cell(size(Value),Model.N);
    else
        BlankValue = zeros([size(Value) Model.N]);
    end
    JackKnife = setfield(JackKnife,FieldNames{i},BlankValue);
end

% perform the jack-knife leave one out procedure
for i = 1:Model.N
    % leave one out
    Include = 1:Model.N;
    Include(i) = 0;
    tempModel.data = Model.data(find(Include),:);
    tempResults = WIPsubfnFitModel(tempModel);
    % store all the parameter estimates in the JackKnife structure
    JackKnife.beta(:,:,i) = tempResults.beta;
    JackKnife.B(:,:,i) = tempResults.B;
    JackKnife.Paths{i} = tempResults.Paths;
end

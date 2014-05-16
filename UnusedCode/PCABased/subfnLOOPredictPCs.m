function predictedY = subfnLOOPredictPCs(data,beta,ModelNum,predictedM,NselPCs,SubjectIndex)

% predictedM: the predicted value of the combination of SSFs for this ONE
% person
% NselPCs: The number of PCs used to create this pattern
% data: the data for this All subjects
% beta: the regression coefficients derived from the rest of the sample
NCOV = size(data.COV,2);
switch ModelNum
    case '4'
        % standard linear regression
        % Expected: 
        % M: [NSub x one or more SSFs]
        % X: [NSub x 1]
        % Y: [NSub x 1]
        % COV: [NSub x 0 or more]
        if NCOV
            predictedY = 1*beta(1) + predictedM*beta(2) + [data.X(SubjectIndex,:) data.COV(SubjectIndex,:)]*beta(3:end);
            
        else
            predictedY = 1*beta(1) + predictedM*beta(2) + [data.X(SubjectIndex,:)]*beta(3:end);
        end
end
        
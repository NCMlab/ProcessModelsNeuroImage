function outAllData = subfnMissingValues(inAllData,MissingDataValue,MissingMethod)
%% Missing inAllData
if nargin == 2
    MissingMethod = 'CasewiseRemoval';
end
outAllData = inAllData;
switch MissingMethod
    case 'MeanReplace'
        % mean replacement
        for i = 1:size(inAllData.M,2)
                outAllData.M(:,i,:) = subfnMeanReplace(squeeze(inAllData.M(:,i,:)),MissingDataValue);
        end
        for i = 1:size(inAllData.X,2)
            outAllData.X(:,i) = subfnMeanReplace(inAllData.X(:,i),MissingDataValue);
        end
        for i = 1:size(inAllData.Y,2)
            outAllData.Y(:,i) = subfnMeanReplace(inAllData.Y(:,i),MissingDataValue);
        end
        if ~isempty(inAllData.COV)
            outAllData.COV = subfnMeanReplace(inAllData.COV,MissingDataValue);
        end
        if ~isempty(inAllData.V)
            outAllData.V = subfnMeanReplace(inAllData.V,MissingDataValue);
        end
        if ~isempty(inAllData.W)
            outAllData.W = subfnMeanReplace(inAllData.W,MissingDataValue);
        end
        if ~isempty(inAllData.Q)
            outAllData.Q = subfnMeanReplace(inAllData.Q,MissingDataValue);
        end
        if ~isempty(inAllData.R)
            outAllData.R = subfnMeanReplace(inAllData.R,MissingDataValue);
        end
       
    case 'CasewiseRemoval'
        NSub = size(inAllData.X,1);
        missF = [];
        for i = 1:size(inAllData.X,2)
            missF = subfnFindMissing(inAllData.X,MissingDataValue,missF);
        end
        for i = 1:size(inAllData.Y,2)
            missF = subfnFindMissing(inAllData.Y,MissingDataValue,missF);
        end
        for i = 1:size(inAllData.M,2)
            for j = 1:size(inAllData.M,3)
                missF = subfnFindMissing(squeeze(inAllData.M(:,i,j)),MissingDataValue,missF);
            end
        end
        if ~isempty(inAllData.COV)
            for i = 1:size(inAllData.COV,2)
               missF = subfnFindMissing(squeeze(inAllData.COV(:,i)),MissingDataValue,missF);
            end
        end
        if ~isempty(inAllData.V)
            missF = subfnFindMissing(inAllData.V,MissingDataValue);
        end
        if ~isempty(inAllData.W)
            missF = subfnFindMissing(inAllData.W,MissingDataValue);
        end
        if ~isempty(inAllData.Q)
            missF = subfnFindMissing(inAllData.Q,MissingDataValue);
        end
        if ~isempty(inAllData.R)
            missF = subfnFindMissing(inAllData.R,MissingDataValue);
        end
        
        missF = unique(missF);
        foundF = ones(NSub,1);
        foundF(missF) = 0;
        includeF = find(foundF);
        outAllData.X = inAllData.X(includeF,:);
        outAllData.M = inAllData.M(includeF,:,:);
        outAllData.Y = inAllData.Y(includeF,:);
        if ~isempty(inAllData.COV)
            outAllData.COV = inAllData.COV(includeF,:);
        end
        if ~isempty(inAllData.V)
            outAllData.V = inAllData.V(includeF,1);
        end
        if ~isempty(inAllData.W)
            outAllData.W = inAllData.W(includeF,1);
        end
        if ~isempty(inAllData.Q)
            outAllData.Q = inAllData.Q(includeF,1);
        end
        if ~isempty(inAllData.R)
            outAllData.R = inAllData.R(includeF,1);
        end
        if ~isempty(inAllData.STRAT)
            outAllData.STRAT = inAllData.STRAT(includeF,1);
        end
end






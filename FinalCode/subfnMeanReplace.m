function outData = subfnMeanReplace(inData,MissingDataValue)
outData = inData;
for i = 1:size(inData,2)
    missF = find(inData(:,i) == MissingDataValue);
    foundF = find(inData(:,i) ~= MissingDataValue);
    % mean replacement
    outData(missF,i) = mean(inData(foundF,i));
end
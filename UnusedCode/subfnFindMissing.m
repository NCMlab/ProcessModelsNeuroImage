function missF = subfnFindMissing(inData,MissingDataValue,missF)
% append any new missing values to the input list
for i = 1:size(inData,2)
    missF = [missF; find(inData(:,i) == MissingDataValue)];
end
function outStruct = subfnSetParameters(name, statsStruct, index)

str = [name ' = {};'];
eval(str);
str = [name ' = setfield(' name ',''beta'', ' num2str(statsStruct.beta(index)) ');'];
eval(str)
str = [name ' = setfield(' name ',''se'', ' num2str(statsStruct.tstat.se(index)) ');'];
eval(str)
str = [name ' = setfield(' name ',''t'', ' num2str(statsStruct.tstat.t(index)) ');'];
eval(str)
str = [name ' = setfield(' name ',''p'', ' num2str(statsStruct.tstat.pval(index)) ');'];
eval(str)
str = ['outStruct = ' name ';'];
eval(str)
% Add standardized beta
str = [name ' = setfield(' name ',''B'', ' num2str(statsStruct.beta(index)*sqrt(statsStruct.covb(index,index))/sqrt(statsStruct.covb(1,1))) ');'];
eval(str)

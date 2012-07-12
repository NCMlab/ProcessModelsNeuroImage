function Model = subfnSetModelParameters(statsStruct)

Model = {};
str = ['Model = setfield(Model,''rsquare'', ' num2str(statsStruct.rsquare) ');'];
eval(str)
str = ['Model = setfield(Model,''adjrsquare'', ' num2str(statsStruct.adjrsquare) ');'];
eval(str)
str = ['Model = setfield(Model,''F'', ' num2str(statsStruct.fstat.f) ');'];
eval(str)
str = ['Model = setfield(Model,''df1'', ' num2str(statsStruct.fstat.dfr) ');'];
eval(str)
str = ['Model = setfield(Model,''df2'', ' num2str(statsStruct.fstat.dfe) ');'];
eval(str)
str = ['Model = setfield(Model,''p'', ' num2str(statsStruct.fstat.pval) ');'];
eval(str)
str = ['Model = setfield(Model,''AIC'', ' num2str(statsStruct.AIC) ');'];
eval(str)



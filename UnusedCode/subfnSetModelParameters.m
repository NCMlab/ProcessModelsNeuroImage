function Model = subfnSetModelParameters(statsStruct)
if isfield(statsStruct,'modelLL')
    % This was logistic regression model
    Model = {};
    str = ['Model = setfield(Model,''modelLL'', ' num2str(statsStruct.modelLL) ');'];
    eval(str)
    str = ['Model = setfield(Model,''modelLLdf'', ' num2str(statsStruct.modelLLdf) ');'];
    eval(str)
    str = ['Model = setfield(Model,''modelLLp'', ' num2str(statsStruct.modelLLp) ');'];
    eval(str)
    str = ['Model = setfield(Model,''m2LL'', ' num2str(statsStruct.m2LL) ');'];
    eval(str)
    str = ['Model = setfield(Model,''McFadden'', ' num2str(statsStruct.McFadden) ');'];
    eval(str)
    str = ['Model = setfield(Model,''CoxSnell'', ' num2str(statsStruct.CoxSnell) ');'];
    eval(str)
    str = ['Model = setfield(Model,''Nagelkrk'', ' num2str(statsStruct.Nagelkrk) ');'];
    eval(str)
    str = ['Model = setfield(Model,''propCorClass'', ' num2str(statsStruct.propCorClass) ');'];
    eval(str)
else
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
    str = ['Model = setfield(Model,''LOOCV'', ' num2str(statsStruct.CV) ');'];
    eval(str)
end
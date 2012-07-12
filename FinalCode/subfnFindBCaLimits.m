function [Alpha1 Alpha2] = subfnFindBCaLimits(bstat,pointEst,alpha,data)
[N Nmed] = size(data.M);
NParameters = size(bstat,3);
NCov = size(data.COV,2);
nboot = length(bstat);


% perform jackknife
theta = zeros(N,Nmed,NParameters);
EMPTYtemp = {};
EMPTYtemp.Y = [];
EMPTYtemp.M = [];
EMPTYtemp.X = [];
EMPTYtemp.V = [];
EMPTYtemp.W = [];
EMPTYtemp.Q = [];
EMPTYtemp.R = [];
EMPTYtemp.COV = [];
EMPTYtemp.Thresholds = data.Thresholds;
EMPTYtemp.ModelNum = data.ModelNum;
if isfield(data,'ProbeMod')
    EMPTYtemp.ProbeMod = data.ProbeMod;
end
Alpha1 = zeros(Nmed,NParameters);
Alpha2 = zeros(Nmed,NParameters);
for i = 1:N
    Values = ones(N,1);
    Values(i) = 0;
    temp = EMPTYtemp;
    % resample the data, NOTE this needs to be expanded for 
    temp.Y = data.Y(logical(Values));
    temp.X = data.X(logical(Values));
    if ~isempty(data.V);temp.V = data.V(logical(Values));end
    if ~isempty(data.W);temp.W = data.W(logical(Values));end
    if ~isempty(data.Q);temp.Q = data.Q(logical(Values));end
    if ~isempty(data.R);temp.R = data.R(logical(Values));end
    for j = 1:Nmed
        temp.M(:,j) = data.M(logical(Values),j);
    end
    for j = 1:NCov
        temp.COV(:,j) = data.COV(logical(Values),j);
    end
    tParam = subfnProcessModelFit(temp,0);
    theta(i,:,:) = tParam.values;
    
end

% 
for k = 1:NParameters
    for j = 1:Nmed
        zh0 = norminv(length(find(bstat(:,j,k) < pointEst(j,k)))/nboot);
        zA = norminv(alpha/2);
        z1mA = norminv(1 - alpha/2);
        ThetaDiff = (sum(theta(:,j))/N) - theta(:,j);
        acc = (sum(ThetaDiff.^3))/(6*(sum(ThetaDiff.^2))^(3/2));
        
        Alpha1(j,k) = normcdf(zh0 + (zh0+zA)/(1 - acc*(zh0 + zA)));
        Alpha2(j,k) = normcdf(zh0 + (zh0+z1mA)/(1 - acc*(zh0 + z1mA)));
    end
end

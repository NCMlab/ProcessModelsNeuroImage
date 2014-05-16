function [Alpha1 Alpha2 p Z] = subfnFindBCaLimits(bstat,pointEst,alpha,data)
%
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
p = zeros(Nmed,NParameters);
Z = zeros(Nmed,NParameters);
% Here is the jack-knife step
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
testValue = 0;
for k = 1:NParameters
    for j = 1:Nmed
        zh0 = norminv(length(find(bstat(:,j,k) < pointEst(j,k)))/nboot);
        zA = norminv(alpha/2);
        z1mA = norminv(1 - alpha/2);
        ThetaDiff = (sum(theta(:,j))/N) - theta(:,j);
        % Calculate the acceleration factor
        acc = (sum(ThetaDiff.^3))/(6*(sum(ThetaDiff.^2))^(3/2));
        Alpha1(j,k) = normcdf(zh0 + (zh0+zA)/(1 - acc*(zh0 + zA)));
        Alpha2(j,k) = normcdf(zh0 + (zh0+z1mA)/(1 - acc*(zh0 + z1mA)));
        % Find percentile of the distribution below the null value
        PCTlower = sum(bstat(:,j,k) < testValue)./nboot;
        PCTupper = sum(bstat(:,j,k) > testValue)./nboot;
        IsLower = PCTlower < PCTupper;
        
        
        % Check to make sure the tails are not zero
        PCTlower = max(PCTlower,1./nboot);
        PCTupper = min(PCTupper, 1-1./nboot);
        
        PCTupper = max(PCTupper,1/nboot);
        PCTlower = min(PCTlower,1-1/nboot);
        
        % Z-score
        ZPCTlower = norminv(PCTlower) - zh0;
        ZPCTupper = norminv(PCTupper) - zh0;
        % adjust for acceleration
        Zlower = (ZPCTlower.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTlower);
        Zupper = (ZPCTupper.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTupper);
        % which tail are we in?
        
        ZTemp = Zupper;
        ZTemp(IsLower) = -Zlower(IsLower);
        Z(j,k) = ZTemp;
        
        pTemp = normcdf(ZTemp);
        pTemp = [pTemp; 1 - pTemp];
        p(j,k) = 2.*min(pTemp);
        
        
    end
end

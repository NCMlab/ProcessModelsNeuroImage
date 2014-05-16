function [bstat, k2stat] = subfnBootStrpv2(data)
% run once to find our how many values get bootstrapped
[ParameterToBS] = subfnProcessModelFit(data,0);
% find the size of the different variables
[Nmed, NParameters] = size(ParameterToBS.values);

N = length(data.Y);
NCov = size(data.COV,2);
% initialize the output values
Nboot = size(data.Resamples,2);
bstat = zeros(Nboot,Nmed,NParameters);
k2stat = zeros(Nboot,1);
% find the number of groups if stratified resampling will be done
% start the bootstrap loop using parallel processing
Resamples = data.Resamples;
parfor i = 1:Nboot
    Samp = Resamples(:,i);
    temp = data;
    % resample the data, NOTE this needs to be expanded for 
    temp.Y = temp.Y(Samp);
    temp.X = temp.X(Samp);
    if ~isempty(temp.V);temp.V = temp.V(Samp);end
    if ~isempty(temp.W);temp.W = temp.W(Samp);end
    if ~isempty(temp.Q);temp.Q = temp.Q(Samp);end
    if ~isempty(temp.R);temp.R = temp.R(Samp);end
    for j = 1:Nmed
        temp.M(:,j) = temp.M(Samp,j);
    end
    for j = 1:NCov
        temp.COV(:,j) = temp.COV(Samp,j);
    end
    tParam = subfnProcessModelFit(temp,0);
    bstat(i,:,:) = tParam.values;
    k2stat(i) = tParam.k2;
end


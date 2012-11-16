function [bstat, k2stat] = subfnBootStrp(data,Nboot)
% run once to find our how many values get bootstrapped
[ParameterToBS] = subfnProcessModelFit(data,0);
% find the size of the different variables
[Nmed, NParameters] = size(ParameterToBS.values);

N = length(data.Y);
NCov = size(data.COV,2);
% initialize the output values
bstat = zeros(Nboot,Nmed,NParameters);
k2stat = zeros(Nboot,1);
% find the number of groups if stratified resampling will be done
if ~isempty(data.STRAT)
    Gr1 = find(data.STRAT == 0);
    NGr1 = length(Gr1);
    Gr2 = find(data.STRAT == 1);
    NGr2 = length(Gr2);
else
    NGr1 = [];
    NGr2 = [];
end
% start the bootstrap loop using parallel processing
for i = 1:Nboot
    % this is needed for the parallel processing to work
    temp = data;
    % create the resamples
    if isempty(temp.STRAT)
        Samp =  floor(N*rand(N,1))+1;
    else
        Samp1 = floor(NGr1*rand(NGr1,1))+1;
        Samp2 = floor(NGr2*rand(NGr2,1))+1+NGr1;
        Samp = [Samp1; Samp2];
    end
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


function outdata = subfnReturnBootStrapData(indata, samples)
if nargin == 1
    N = length(indata.Y);
    % find the number of groups if stratified resampling will be done
    if ~isempty(indata.STRAT)
        Gr1 = find(indata.STRAT == 0);
        NGr1 = length(Gr1);
        Gr2 = find(indata.STRAT == 1);
        NGr2 = length(Gr2);
    else
        NGr1 = [];
        NGr2 = [];
    end
    outdata = indata;
    if isempty(indata.STRAT)
        Samp =  floor(N*rand(N,1))+1;
    else
        Samp1 = floor(NGr1*rand(NGr1,1))+1;
        Samp2 = floor(NGr2*rand(NGr2,1))+1+NGr1;
        Samp = [Samp1; Samp2];
    end
end
if isfield(indata,'COV')
    NCov = size(indata.COV,2);
else
    NCov = 0;
end
Nmed = size(indata.M,2);
% resample the data, NOTE this needs to be expanded for
outdata.Y = indata.Y(Samp);
outdata.X = outdata.X(Samp);
if ~isempty(outdata.V);outdata.V = outdata.V(Samp);end
%if ~isempty(outdata.W);outdata.W = outdata.W(Samp);end
%if ~isempty(outdata.Q);outdata.Q = outdata.Q(Samp);end
%if ~isempty(outdata.R);outdata.R = outdata.R(Samp);end
for j = 1:Nmed
    outdata.M(:,j) = outdata.M(Samp,j);
end
for j = 1:NCov
    outdata.COV(:,j) = outdata.COV(Samp,j);
end
function Samp = WIPBootStrap(STRAT,N)
% For each cluster chunk it is creating the EXACT same set of bootstrap
% resamples ... this is bad. The following line attempts to avoid that by
% resetting the random seed.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
if ~isempty(STRAT)
    Gr1 = find(STRAT == 0);
    NGr1 = length(Gr1);
    Gr2 = find(STRAT == 1);
    NGr2 = length(Gr2);
else
    NGr1 = [];
    NGr2 = [];
end

% create the resamples
if isempty(STRAT)
    %N = size(data.data,1);
    Samp =  floor(N*rand(N,1))+1;
else
    Samp1 = floor(NGr1*rand(NGr1,1))+1;
    Samp2 = floor(NGr2*rand(NGr2,1))+1+NGr1;
    Samp = [Samp1; Samp2];
end
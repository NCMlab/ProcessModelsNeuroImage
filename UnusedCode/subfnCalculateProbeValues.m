function output = subfnCalculateProbeValues(input)
% Calculate the probe values at the following percentiles:
% 10, 25, 50, 75 90 
% and at plus or minus one standard deviation
mIn = mean(input);
stdIn = std(input);

output = [0 sort([prctile(input,[10 25 50 75 90]) mIn-stdIn mIn+stdIn])];

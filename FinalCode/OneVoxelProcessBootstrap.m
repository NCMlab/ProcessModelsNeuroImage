function [Results, BootStrap] = OneVoxelProcessBootstrap(Model)
% reset the random number generator. This is EXTREMEMLY important when
% deploying these analyses to a cluster. Without this resetting then each
% node of the cluster CAN choose the exact same random numbers.

rng('shuffle','multFibonacci');
%rng(seed,'twister')
t = tic;
% Fit the model
% The indirect paths (the Paths cell) is set to be an array of cells. This
% is rather a pain in the neck for everything else. 
Results = FitProcessModel(Model);
% perform bootstrap
FieldNames = {'beta' 'B' };%'Paths'};
if Model.Nboot > 0
    Model.STRAT = [];
    % The bootstrap propcedure returns the bootstrap distributions from the
    % field names stated above. These same effects are then run through
    % the Jack-Knife procedure.
    try 
       % Perform the jack-knife step
       
        JackKnife = JackKnifeFunction(Model,Results,FieldNames);    
        BootStrap = BootStrapFunction(Model,Model.Nboot,FieldNames);
        
        % Calculate BCaCI for beta values
        [Z, p] = CalculateBCaCIsimple(Results.beta,BootStrap.beta, JackKnife.beta, Model.Thresholds(1));
        % put these values into the output structure
        Results.BCaCI.Z = Z;
        Results.BCaCI.p = p;
        % Calculate BCaCI for Path Values
        % Cycle over probes in bootstrap values
        Nprobe = length(Results.ProbeValues);
        BCaCIpathZ = zeros(Nprobe,1);
        BCaCIpathp = zeros(Nprobe,1);
        for j = 1:Nprobe
            % extract bootstrap values for each probe vcalue
            tempBootstrap = zeros(Model.Nboot,1);
            for k = 1:Model.Nboot
                tempBootstrap(k) = BootStrap.Paths{k}{1}{1}(j);
            end
            % extract JackKnife values for each probe vcalue
            tempJackKnife = zeros(Model.Nsub,1);
            for k = 1:Model.Nsub
                tempJackKnife(k) = JackKnife.Paths{k}{1}{1}(j);
            end
            [BCaCIpathZ(j), BCaCIpathp(j)] = CalculateBCaCIsimple(Results.Paths{1}{1}(j),tempBootstrap,tempJackKnife, Model.Thresholds(1));
        end  
        Results.BCaCIpathZ = BCaCIpathZ;
        Results.BCaCIpathp = BCaCIpathp;
        
        
        %[Alpha1 Alpha2 Z p] = CalculateBCaLimitsOneValue(squeeze(JackKnife.beta(2,2,:)),Results.beta(2,2), squeeze(BootStrap.beta(2,2,:)),0.05);
        
        
%        BoLBBcACI = BagOfBootStrapFunction(Model,100,FieldNames,JackKnife);
        % Calculate the BCaci values for each parameter
 %%%%%%       Results.BCaCI = CreateBCaCI(Results,BootStrap,JackKnife,Model.Thresholds);
%        Results.BoLBBCaCI = BoLBBcACI;
    catch me
%        error('error!');
        Results = [];
    end
end
%toc(t);

%% TO DO: unflatten the path BCaCI results

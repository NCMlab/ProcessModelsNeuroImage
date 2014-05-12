function Results = WIPsubfnProcessData(Model)
% reset the random number generator. This is EXTREMEMLY important when
% deploying these analyses to a cluster. Without this resetting then each
% node of the cluster CAN choose the exact same random numbers.
%RandStream.setDefaultStream(RandStream('mt19937ar','Seed',sum(100*clock)));

% Fit the model
% The indirect paths (the Paths cell) is set to be an array of cells. This
% is rather a pain in the neck for everything else. 
Results = WIPsubfnFitModel(Model);
% perform bootstrap
FieldNames = {'beta' 'B' 'Paths'};
if Model.Nboot > 0
    Model.STRAT = [];
    % The bootstrap propcedure returns the bootstrap distributions from the
    % field names stated above. These same effects are then run through
    % the Jack-Knife procedure.
    BootStrap = WIPsubfnBootStrap(Model,Model.Nboot,FieldNames);

    % Perform the jack-knife step
    JackKnife = WIPJackKnife(Model,Results,FieldNames);
    
    % Calculate the BCaci values for each parameter
    Results.BCaCI = WIPsubfnCreateBCaCI(Results,BootStrap,JackKnife,Model.Thresh);
  
end


%% TO DO: unflatten the path BCaCI results

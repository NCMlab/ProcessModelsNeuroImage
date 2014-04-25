function Results = WIPsubfnProcessData(Model)
RandStream.setDefaultStream(RandStream('mt19937ar','Seed',sum(100*clock)));
% Fit the model
Results = WIPsubfnFitModel(Model);
% perform bootstrap
if Model.Nboot > 0
    Model.STRAT = [];
    BootStrap = WIPsubfnBootStrap(Model,Model.Nboot);
    % Calculate BCaci for each beta and the path
    
    % Perform the jack-knife step

    JackKnife = WIPJackKnife(Model,Results);
    
    
    % Calculate the BCaci values for each parameter
    % beta
    BootStrapData = BootStrap.beta;
    JackKnifeData = JackKnife.beta;
    PointEstimate = Results.beta;
    [Alpha1 Alpha2 Z p] = WIPsubfnCalculateBCaci(JackKnifeData,PointEstimate, BootStrapData,alpha);
    % To Do: [Sbstat(ceil(Alpha1(j,k)*nboot),j,k) Sbstat(ceil(Alpha2(j,k)*nboot),j,k)]
    % indirect effects
    [m n] = size(BootStrap.Paths{1});
    BootStrapData = zeros(m,n,Model.Nboot);
    for i = 1:Model.Nboot
        BootStrapData(:,:,i) = BootStrap.Paths{i}{:};
    end
    JackKnifeData = zeros(m,n,Model.N);
    for i = 1:Model.N
        JackKnifeData(:,:,i) = JackKnife.Paths{i}{:};
    end
    PointEstimate = Results.Paths{:};
    [Alpha1 Alpha2 Z p] = WIPsubfnCalculateBCaci(JackKnifeData,PointEstimate, BootStrapData,alpha);

end
% perform permutation
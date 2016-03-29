function [OutPath, OutSE, probeValues] = subfnCalculatePathSE(Results, data)
[N M] = size(data.data);
Direct = data.Direct;
Inter = data.Inter;
Paths = data.Paths;
probeValues = [];
% cycle over the different paths
for j = 1:size(Paths,3)
    % cycle over the steps in the path
   
    NSteps = max(max(Paths(:,:,j)));
    % Create a cell array where each cell is the parameter(s) for that
    % step. The cells can be arrays if the path is modulated.
    PathsParameters = cell(NSteps,1);
    PathsStandardErrors = cell(NSteps,1);
    for i = 1:max(max(Paths(:,:,j)))
        
        % Which interaction parameter is IN the path?
        % Which is NOT
        
        index = find(Paths(:,:,j)==i);
        [Row Col] = ind2sub(size(Paths),index);
        % is there an interaction on this step?
        tempInter = data.Inter(:,Col);
        if isempty(find(tempInter))
            % No interaction for this model
            PathsParameters{i} = Results.beta(Row+1,Col);
            PathsStandardErrors{i} = (Results.covb(Row+1,Row+1,Col));
        else
            % YES, there is an interaction in this model?
            F = find(tempInter) + 1;
            % Is the interaction effect PART OF THIS Path?
            if ~isempty(find(F == Row + 1)) % CHANGED Col to Row + 1
                % YES it is
                % which variable do you probe?
                F(find(F == Row + 1)) = 0; % MADE SAME CHANGE HERE
                % The following are in "parameter space"
                ParameterToProbe = F(find(F));
                DirectEffectParameter = Row + 1;
                % Find the data variable to probe
                Fmod = F(find(F)) - 1;
                % do not probe a dichotomous variable
                
                Moderators = data.data(:,Fmod);
               % find the probe values,this even works for higher
               % order interactions, I think.
                if length(unique(Moderators)) == 2
                    probeMod = unique(Moderators);
                    probeValues = probeMod;
                else
                    %probeMod = prctile(Moderators,[10:10:90]);
                    probeMod = prctile(Moderators,[5:5:95]);
                    probeValues = [zeros(size(probeMod,1),1) probeMod];
                    %  probeValues = probeValues(:,1)*probeValues(:,2)';
                end
                % What is the effect for this step of the path?
                % It is the parameter for the direct effect PLUS the
                % interaction effect multiplied by the moderator at the
                % probed values.
                PathsParameters{i} = Results.beta(DirectEffectParameter,Col) + Results.beta(M+1+1,Col).*probeValues;
                % Get the standard error of the direct component
                % Standard error of the interaction term
                PathsStandardErrors{i} = Results.covb(DirectEffectParameter,DirectEffectParameter,Col) + ...% good
                2.*Results.covb(M+2,DirectEffectParameter,Col).*probeValues + ... % I am not sure if the use of M will generalize
                Results.covb(M+2,M+2,Col).*probeValues.^2;
            else
                PathsParameters{i} = Results.beta(Row+1,Col);
                PathsStandardErrors{i} = (Results.covb(Row+1,Row+1,Col));
            end
        end
    end
end
     
Path = 1;
for i = 1:length(PathsParameters)
    Path = PathsParameters{i}.*Path;
end

OutPath = Path;
OutSE = 0;
for i = 1:NSteps
    AllSteps = 1:NSteps;
    CurrentStep = i;
    AllSteps(i) = 0;
    OtherSteps = AllSteps(find(AllSteps));
    ThisStepSE = 1;
    for j = 1:length(OtherSteps)
        ThisStepSE = ThisStepSE.*PathsParameters{OtherSteps(j)}.^2;
    end
    ThisStepSE = ThisStepSE.*PathsStandardErrors{CurrentStep};
    OutSE = OutSE + ThisStepSE;
end
OutSE = sqrt(OutSE);
    
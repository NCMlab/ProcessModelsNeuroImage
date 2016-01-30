function subfnCalculatePathSE(Results, Direct, Inter, Paths)

Direct = data.Direct;
Inter = data.Inter;
Paths = data.Paths;
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
            PathsParameters{i} = Results.beta(Row,Col);
            PathsStandardErrors{i} = sqrt(Results.covb(Row,Row,Col));
        else
            % YES, there is an interaction in this model?
            F = find(tempInter) + 1;
            % Is the interaction effect PART OF THIS Path?
            if ~isempty(find(F == Col))
                % YES it is
                % which variable do you probe?
            
            % find the probe values,this even works for higher
            % order interactions, I think.
            
          
                F(find(F == Col)) = 0;
                % do not probe a
                Fmod = F(find(F)) - 1;
                Moderators = data.data(:,Fmod);
                if length(unique(Moderators)) == 2
                    probeMod = unique(Moderators);
                    probeValues = probeMod;
                else
                    %probeMod = prctile(Moderators,[10:10:90]);
                    probeMod = prctile(Moderators,[5:5:95]);
                    probeValues = [zeros(size(probeMod,1),1) probeMod];
                    %  probeValues = probeValues(:,1)*probeValues(:,2)';
                end
                InteractionComponent = InteractionComponent + Results.beta(M+1+1,Col).*probeValues;
                ResultPath = ResultPath.*InteractionComponent;
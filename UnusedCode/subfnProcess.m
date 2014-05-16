function Parameters = subfnProcess(temp)

% it initially does not probe the moderator

[pointEst Parameters{1}] = subfnProcessModelFit(temp,1);
% here we check to see if the interaction is significant.
% if so then probe the mod for all bootstraps
if pointEst.ProbeMod
    temp.ProbeMod = 1;
    % since the interaction is significant then we want to probe
    % the moderator by re-running the regression
    [pointEst Parameters{1}] = subfnProcessModelFit(temp,1);
end

Nsub = size(temp.Y,1);
% have the option of turning the boot strapping off.
% turn off boot strapping on the interaction effect if the
% Johnson-Neyman shows no range of significant interactions.
%
if temp.Nboot %& Parameters{i}.JohnsonNeyman ~= -99
    % calculate the boot strap estimates
    [bstat k2stat] = subfnBootStrp(temp,temp.Nboot);
    %[bstat k2stat] = subfnBootStrpv2(temp);
    if k2stat(1) ~= 0
        Sk2stat = sort(k2stat);
        k2 = {};
        k2.pointEst = pointEst.k2;
        PERci = cell(length(temp.Thresholds));
        
        for k = 1:length(temp.Thresholds)
            t = num2str(temp.Thresholds(k));
            PERci{k} = setfield(PERci{k},['alpha' t(3:end)],...
                [Sk2stat(ceil(temp.Thresholds(1)/2*temp.Nboot)) Sk2stat(ceil((1-temp.Thresholds(1)/2)*temp.Nboot))]);
        end
        k2.PERci = PERci;
        Parameters{1}.k2 = k2;
    end
    Nsub = size(temp.X,1);
    [nboot Nmed NParameters] = size(bstat);
    % Find the number of non-zero bstat values
    if length(find(squeeze(bstat(1,1,:)))) == 1
        NParameters = 1;
    end
    
    % calculate the Confidence intervals on the parameters that
    % need to be bootstrapped.
    [BCaci PERci p Z] = subfnFindConfidenceIntervals(temp,bstat,pointEst);
    % pValue = length(find(pointEst.values>bstat))/temp.Nboot;
     
%     pValue = zeros(length(pointEst.values),1);
%     for kk = 1:length(pointEst.values)
%         pValue(kk) = length(find(pointEst.values(kk)>bstat(:,1,kk)))/temp.Nboot;
%     end
%    pointEst.pValue = pValue';
    % Then fill in the appropriate Parameters with the confidence
    % intervals.
    %            str = [pointEst.names ' = {};'];
    %            eval(str)
    
    for j = 1:Nmed
        temp2 = {};
        for k = 1:NParameters
            %                     str = sprintf('temp2 = Parameters{i}.%s{j};',pointEst.names);
            %                     eval(str);
            str = sprintf('temp2=setfield(temp2,''BCaci'',BCaci{%d,%d});',j,k);
            eval(str);
            str = sprintf('temp2=setfield(temp2,''PERci'',PERci{%d,%d});',j,k);
            eval(str);
            str = sprintf('temp2=setfield(temp2,''bootSE'',bootSE(bstat(:,%d,%d),Nsub));',j,k);
            eval(str);
            str = sprintf('temp2=setfield(temp2,''pointEst'',pointEst.values(j,k));');
            eval(str);
             str = sprintf('temp2=setfield(temp2,''p'',p(j,k));');
             eval(str);
             str = sprintf('temp2=setfield(temp2,''Z'',Z(j,k));');
             eval(str);
%             str = sprintf('temp2=setfield(temp2,''pValue'',pointEst.pValue);');
%             eval(str);
            if isfield(pointEst,'probeValues')
                str = sprintf('temp2=setfield(temp2,''probeValue'',pointEst.probeValues(k));');
                eval(str);
            end
            str = sprintf('%s{%d}=temp2;',pointEst.names(j,:),k);
            eval(str);
        end
        str = sprintf('Parameters{1}.%s = %s;',pointEst.names(j,:),pointEst.names(j,:));
        eval(str);
    end
    Parameters{1}.Nboot = temp.Nboot;
    Parameters{1}.Thresholds = temp.Thresholds;
else
    for j = 1:size(pointEst,1);
        temp2 = {};
        for k = 1:size(pointEst,2);
            str = sprintf('temp2=setfield(temp2,''BCaci'',[]);');
            eval(str);
            str = sprintf('temp2=setfield(temp2,''PERci'',[]);');
            eval(str);
            str = sprintf('temp2=setfield(temp2,''bootSE'',[]);');
            eval(str);
            str = sprintf('%s{%d,%d}=temp2;',pointEst.names,k,j);
            eval(str);
        end
    end
end
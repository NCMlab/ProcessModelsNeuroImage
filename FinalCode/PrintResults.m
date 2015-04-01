function PrintResults(Model,Results)
%%
fid = 1;
fprintf(fid,'-------------------------------------------\n');
NVar = size(Model.data,2);
% print each model
for i = 1:NVar
    if sum(Model.Direct(:,i))
        % print dependent variable
        fprintf(fid,'\nDependent variable: %s\n',Model.Names{i});
        % print independent variables
        Rows = find(Model.Direct(:,i));
        %fprintf(fid,'%20s\tUn-Coef\t\tSt-Coef\t\t\n','');
        fprintf(fid,'%20s%10s%10s%10s%10s%10s\n','Model','beta','stderr','B','t','p');
        fprintf(fid,'%20s','constant');
        fprintf(fid,'%10.3f%10.3f%10s%10.3f%10.3f\n',...
            Results.beta(1,i),...
            Results.se(1,i),...
            '--',...
            Results.t(1,i),...
            Results.p(1,i));
        for j = 1:length(Rows)
            fprintf(fid,'%20s',Model.Names{Rows(j)});
            fprintf(fid,'%10.3f%10.3f%10.3f%10.3f%10.3f\n',...
                Results.beta(Rows(j)+1,i),...
                Results.se(Rows(j)+1,i),...
                Results.B(Rows(j)+1,i),...
                Results.t(Rows(j)+1,i),...
                Results.p(Rows(j)+1,i));
        end
    end
    % Are there any interactions for this variable?
    if sum(Model.Inter(:,i))
        % How many interactions are there in this model?
        NInter = size(Results.beta,1) - (Model.Nvar+1);
        for k = 1:NInter
            % Create the name of the interaction
            InterRows = find(Model.Inter(:,i));
            InterName = '';
            for j = 1:length(InterRows)
                InterName = sprintf('%s%s',InterName,Model.Names{InterRows(j)});
                if j < length(InterRows)
                    InterName = sprintf('%s_x_',InterName);
                end
            end
            % Print out the parameter estimates for the interaction
            fprintf(fid,'%20s',InterName);
            fprintf(fid,'%10.3f%10.3f%10.3f%10.3f%10.3f\n',...
                Results.beta(Model.Nvar+1+k,i),...
                Results.se(Model.Nvar+1+k,i),...
                Results.B(Model.Nvar+1+k,i),...
                Results.t(Model.Nvar+1+k,i),...
                Results.p(Model.Nvar+1+k,i));
        end
    end
end




% Print paths
NPaths = size(Model.Paths,3);
%%

for i = 1:NPaths
    % cycle over multiple paths
    fprintf(fid,'\nPath: %d\n',i);
    NSteps = max(max(Model.Paths(:,:,i)));
    for j = 1:NSteps
        % print out the names for each step in the path
        index = find(Model.Paths(:,:,i)==j);
        [Row Col] = ind2sub(size(Model.Paths),index);
        fprintf(fid,'%20s to %20s\n',Model.Names{Row},Model.Names{Col});
    end
    
    
    for k = 1:length(Model.Thresholds)
        % what is the name of the moderating variable
        % mediator/moderatorValue/effect/BootSE/LLCI/UUCI
%        fprintf(fid,'Indirect Effect size: %0.4f\n',Results.Paths{i});
%        fprintf(fid,'At alpha = %0.4f, BCaCI = [%0.4f, %0.4f]\n',Model.Thresholds(k),Results.BCaCI.Paths(:,:,1,i,k),Results.BCaCI.Paths(:,:,2,i,k));
        fprintf(fid,'alpha = %0.3f\n',Model.Thresholds(k))
        fprintf(fid, 'Number of bootstrap resamples: %d\n',Model.Nboot);
        if length(Results.ProbeValues{i}) > 0
            fprintf(fid,'%10s%10s%10s%10s%10s\n','ProbeVal','EffSize','LLBCa','UUBCa','BCaZ')
            for j = 1:length(Results.ProbeValues{i})
                fprintf(fid,'%10.4f%10.4f%10.3f%10.3f%10.3f\n',Results.ProbeValues{i}(j),Results.Paths{i}(j),...
                    Results.BCaCI.Paths(j,i,1,1,k),...
                    Results.BCaCI.Paths(j,i,2,1,k),...
                    Results.BCaCI.PathsZ(j));
            end
        else
            fprintf(fid,'%10s%10s%10s%10s\n','EffSize','LLBCa','UUBCa','BCaZ');
            fprintf(fid,'%10.4f%10.4f%10.3f%10.3f\n',Results.Paths{i}(1),...
                    Results.BCaCI.Paths(1,1,1,i,k),...
                    Results.BCaCI.Paths(1,1,2,i,k),...
                    Results.BCaCI.PathsZ(i));
        end
    end
end
fprintf(fid,'-------------------------------------------\n');

%% print the correlation matrix of the variables
[r] = corr(Model.data);
fprintf(1,' ==== Correlation Matrix of data ===\n');
fprintf(1,'%10s','');
for i = 1:Model.Nvar
    fprintf(1,'%10s',Model.Names{i});
end
fprintf(1,'\n');
for i = 1:Model.Nvar
    fprintf(1,'%10s',Model.Names{i});
    for j = 1:Model.Nvar
        if j < i
            fprintf(1,'%10.2f',r(i,j));
        else
            fprintf(1,'%10s','-');
        end
    end
    fprintf(1,'\n');
end


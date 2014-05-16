function WIPPrintResults(Model,Results)
%%
fid = 1;
NVar = size(Model.data,2);
% print each model
for i = 1:NVar
    if sum(Model.Direct(:,i))
        % print dependent variable
        fprintf(fid,'\nDependent variable: %s\n',Model.names{i});
        % print independent variables
        Rows = find(Model.Direct(:,i));
        fprintf(fid,'%20s\tUn-Coef\t\tSt-Coef\t\t\n','');
        fprintf(fid,'%20s\tbeta\tstderr\tB\tt\tp\n','Model');
        fprintf(fid,'%20s\t','constant');
        fprintf(fid,'%6.4f\t%6.4f\t\t%6.4f\t%6.4f\n',...
            Results.beta(1,i),...
            Results.se(1,i),...
            Results.t(1,i),...
            Results.p(1,i));
        for j = 1:length(Rows)
            fprintf(fid,'%20s\t',Model.names{j});
            fprintf(fid,'%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
            Results.beta(j+1,i),...
            Results.se(j+1,i),...
            Results.B(j+1,i),...
            Results.t(j+1,i),...
            Results.p(j+1,i));
        end
    end
end
% Print paths
NPaths = size(Model.Paths,3);
for i = 1:NPaths
    fprintf(fid,'\nPath: %d\n',i);
    NSteps = max(max(Model.Paths(:,:,i)));
    for j = 1:NSteps
        index = find(Model.Paths(:,:,i)==j);
        [Row Col] = ind2sub(size(Model.Paths),index);
        fprintf(fid,'%20s to %20s\n',Model.names{Row},Model.names{Col});
    end
    fprintf(fid,'Indirect Effect size: %0.4f\n',Results.Paths{i});
    for k = 1:length(Model.Thresh)
        fprintf(fid,'At alpha = %0.4f, BCaCI = [%0.4f, %0.4f]\n',Model.Thresh(k),Results.BCaCI.Paths(:,:,1,i,k),Results.BCaCI.Paths(:,:,2,i,k));
    end
end

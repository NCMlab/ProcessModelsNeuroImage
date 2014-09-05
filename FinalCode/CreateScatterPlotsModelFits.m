function CreateScatterPlotsModelFits(Model)

NVar = size(Model.data,2);
% Cycle over the steps in the paths
NPaths = size(Model.Paths,3);
%%
Gr1 = find(Model.data(:,1)<0);
Gr2 = find(Model.data(:,1)>0);
for i = 1:NPaths
    f1 = figure;
    f2 = figure;
    clf
    Steps = unique(Model.Paths(:,:,i));
    Steps = Steps(find(Steps));
    NSteps = length(Steps);
    for j = 1:NSteps
        figure(f1)
        [Row Col] = ind2sub(size(Model.Direct),find(Model.Paths(:,:,i) == Steps(j)));
        
        subplot(1,NSteps,j)
        hold on
        scatter(Model.data(Gr1,Row),Model.data(Gr1,Col))
        scatter(Model.data(Gr2,Row),Model.data(Gr2,Col),'r')
        xlabel(Model.Names{Row});
        ylabel(Model.Names{Col});
        % Find the residulaized relationship
        Indep = find(Model.Direct(:,Col));
        Indep(Row) = 0;
        Indep = Indep(find(Indep));
        
        Sind = regstats(Model.data(:,Col),Model.data(:,Indep));
        Sdep = regstats(Model.data(:,Row),Model.data(:,Indep));
        
        figure(f2)
        
        subplot(1,NSteps,j)
        hold on
        scatter(Sdep.r(Gr1),Sind.r(Gr1))
        scatter(Sdep.r(Gr2),Sind.r(Gr2),'r')
        xlabel(Model.Names{Row});
        ylabel(Model.Names{Col});
        S2 = regstats(Sdep.r,Sind.r);
        S2.tstat.t;
        
        S1 = regstats(Model.data(:,Col),Model.data(:,Row));
        S1.tstat.t;
    end
    figure(f1);
    title(sprintf('Path %d',i))
    figure(f2);
    title(sprintf('Path %d',i))
end
    



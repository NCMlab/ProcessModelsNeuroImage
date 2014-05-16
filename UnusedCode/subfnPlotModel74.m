function subfnPlotModel74(Parameters,data)
%%
%resY = data.Y - data.M.*Parameters{1}.Model2.B1.beta - data.X.*Parameters{1}.Model2.A.beta;

[values] = unique(data.X);
NGroups = length(values);
Groups = cell(NGroups,1);
for i = 1:NGroups
    Groups{i} = find(data.X == values(i));
end
Means = zeros(NGroups,2);
Stds = zeros(NGroups,2);
for i = 1:NGroups
    Means(i,1) = mean(data.M(Groups{i}));
    Means(i,2) = mean(data.Y(Groups{i}));
    Stds(i,1) = std(data.M(Groups{i}));
    Stds(i,2) = std(data.Y(Groups{i}));
end
%
colors = {'b' 'r' 'g'};
figure(1)
clf
hold on
for i = 1:NGroups
    plot(data.M(Groups{i}),data.Y(Groups{i}),['o' colors{i}]);
    k = line([Means(i,1)-Stds(i,1) Means(i,1) + Stds(i,1)],[ Means(i,2) Means(i,2)]);
    set(k,'Color',colors{i});
    k = line([Means(i,1) Means(i,1)],[Means(i,2)-Stds(i,2) Means(i,2) + Stds(i,2)]);
    set(k,'Color',colors{i});
    
end
%%
% remove the effects of X on M and Y
S1 = regstats(data.Y,data.X);
S2 = regstats(data.M,data.X);
YmX = S1.r;
MmX = S2.r;

figure(2)
clf
hold on
for i = 1:NGroups
    plot(MmX(Groups{i}),YmX(Groups{i}),['o' colors{i}]);
end




%Parameters{1}
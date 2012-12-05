clear
NSub = 200;
NVox = 1000;
M = randn(NSub,NVox);
NPCs = 12;

Y = randn(NSub,1);

tic
[c1 scores1 latent1] = pca(M);
toc
Age = scores1(:,4) + randn(NSub,1)*0.25 + Y;
Y = scores1(:,[1 3 5 7])*[1 0.4 0.3 0.2]'+randn(NSub,1).*0.25 + Age.*0.25;
corr([Age Y scores1(:,1:7)])
    
fullmodel_beta = regress(Y,[ones(NSub,1) scores1(:,1:NPCs)]);

combo_matrix = boolean_enumeration_f(NPCs);
% how many combos are there?
NCombos = size(combo_matrix,1)

AIC1 = zeros(NCombos,1);
AIC2 = zeros(NCombos,1);
tic
for j = 1:NCombos
    selected_PCs = find(combo_matrix(j,:));
    behav_fit_coef = regress(Y, [ones(NSub,1) scores1(:,selected_PCs) Age]);
    fit1 = [ones(NSub,1) scores1(:,selected_PCs) Age]*behav_fit_coef;
    resid1 = Y - fit1;
    NParam = length(behav_fit_coef);
    AIC1(j) = NSub * log(resid1'*resid1 / NSub )  + ...
        2*(NParam)*(NParam+1)/(NSub-NParam-1) + 2*(NParam);
end
t1=toc
tic
for j = 1:NCombos
    selected_PCs = find(combo_matrix(j,:));
    fit2 = [ones(NSub,1) scores1(:,selected_PCs) Age]*fullmodel_beta([1 selected_PCs+1 end]);
    NParam = length(selected_PCs) + 1;
    resid2 = Y - fit2;
    AIC2(j) = NSub * log(resid2'*resid2 / NSub )  + ...
        2*(NParam)*(NParam+1)/(NSub-NParam-1) + 2*(NParam);
end
t2=toc

t1/t2

corr([AIC1 AIC2])

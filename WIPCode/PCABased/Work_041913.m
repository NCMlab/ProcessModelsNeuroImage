BaseDir = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes';
BaseDir = '/Users/jason/Documents/MyData/ModMedCogRes';
JobsDir = fullfile(BaseDir,'jobs');
DataFile = fullfile(BaseDir,'AllData')
load(DataFile)
data.Thresholds = AnalysisParameters.Thresholds;
if ~exist(JobsDir,'dir')
    % Make the outrput path for cluster jobs
    mkdir(JobsDir)
end

ModelNum = '4';
[NSub NMed NVox] = size(data.M);
if NMed > 1
    error('This doesn''t work with multiple mediators yet!');
end
NCOV = size(data.COV,2);

NPCs = 12;
% Create matrix of all possible cominations of PCs
combo_matrix = boolean_enumeration_f(NPCs);
% how many combos are there?
NCombos = size(combo_matrix,1);
remove_row_means = 1;


% perform PCA
combo_matrix = boolean_enumeration_f(NPCs);
NCombos = size(combo_matrix,1);
[lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(data.M)', remove_row_means);
ssf = squeeze(data.M) * eigenimages_noZeroes;
ssfSubset = ssf(:,1:NPCs);

% Use the set of SSFs as predictors of a vector of ones

AllMedEffects = zeros(NCombos,1);
AIC = zeros(NCombos,3);
Tstats = zeros(N,9);
N = length(data.X);
Yind = find(data.X == 0);
Oind = find(data.X == 1);

for j = 1:NCombos
    selected_PCs = find(combo_matrix(j,:));
    X = ssfSubset(:,selected_PCs);
    % Fits a model of ones
    
    Y = ones(N,1);
    beta = inv(X'*X)*X'*Y;
    fit = X*beta;
    r = Y - fit;
    p = length(beta);
    dfe = N-p;
    AIC(j,1) =  N*log(sum(r.^2)/N) + 2*p*(p+1)/(dfe-1) + 2*p;
    % test to see how the different combos compare across groups
    [H prob CI Stats] = ttest(fit(Yind));
    Tstats(j,1) = Stats.tstat;
    [H prob CI Stats] = ttest(fit(Oind));
    Tstats(j,2) = Stats.tstat;
    [H prob CI Stats] = ttest2(fit(Yind),fit(Oind));
    Tstats(j,3) = Stats.tstat;
    % Fits a YOUNG ONLY model
    Y = ones(N,1);
    Y(Oind) = 0;
    beta = inv(X'*X)*X'*Y;
    fit = X*beta;
    r = Y - fit;
    p = length(beta);
    dfe = N-p;
    AIC(j,2) =  N*log(sum(r.^2)/N) + 2*p*(p+1)/(dfe-1) + 2*p;
    % test to see how the different combos compare across groups
    [H prob CI Stats] = ttest(fit(Yind));
    Tstats(j,4) = Stats.tstat;
    [H prob CI Stats] = ttest(fit(Oind));
    Tstats(j,5) = Stats.tstat;
    [H prob CI Stats] = ttest2(fit(Yind),fit(Oind));
    Tstats(j,6) = Stats.tstat;
    % Fits an OLD ONLY model
    Y = ones(N,1);
    Y(Yind) = 0;
    beta = inv(X'*X)*X'*Y;
    fit = X*beta;
    r = Y - fit;
    p = length(beta);
    dfe = N-p;
    AIC(j,3) =  N*log(sum(r.^2)/N) + 2*p*(p+1)/(dfe-1) + 2*p;
    % test to see how the different combos compare across groups
    [H prob CI Stats] = ttest(fit(Yind));
    Tstats(j,7) = Stats.tstat;
    [H prob CI Stats] = ttest(fit(Oind));
    Tstats(j,8) = Stats.tstat;
    [H prob CI Stats] = ttest2(fit(Yind),fit(Oind));
    Tstats(j,9) = Stats.tstat;
end

% find the best combo for the common pattern for both age groups which is
% used the same
alpha = 0.05;
Tthr = spm_invTcdf(1-alpha,N);
F = find(Tstats(:,1) > Tthr  & Tstats(:,2) > Tthr & Tstats(:,3) > Tthr);
PCs_YngANDOld_DIFF = find(combo_matrix(find(AIC(F,1) == min(AIC(F,1))),:));


F = find(Tstats(:,1) > Tthr & Tstats(:,2) > Tthr & Tstats(:,3) < Tthr);
PCs_YngANDOld_SAME = find(combo_matrix(find(AIC(F,1) == min(AIC(F,1))),:));


F = find(Tstats(:,4) > Tthr & Tstats(:,5) < Tthr & Tstats(:,6) > Tthr);
PCs_YNG = find(combo_matrix(find(AIC(F,2) == min(AIC(F,2))),:));

F = find(Tstats(:,7) < Tthr & Tstats(:,8) > Tthr & Tstats(:,9) > Tthr);
PCs_OLD = find(combo_matrix(find(AIC(F,3) == min(AIC(F,3))),:));

% create the SSFs
% Young AND Old no difference
if ~isempty(PCs_YngANDOld_DIFF)
    X = ssfSubset(:,PCs_YngANDOld_DIFF);
    % Fits a model of ones
    Y = ones(N,1);
    beta = inv(X'*X)*X'*Y;
    SSF_YngANDOld_DIFF = X*beta;
    temp = eigenimages_noZeroes(:, PCs_YngANDOld_DIFF)*beta;
    PC_image_YngANDOld_DIFF = temp / norm(temp);
    clear temp;
    %%%%% Obtain SSFs associated with the normalized best linear behavioral-fit image %%%%%
    PC_image_ssfs_YngANDOld_DIFF = squeeze(data.M) * PC_image_YngANDOld_DIFF;
end
% Young AND Old Difference
if ~isempty(PCs_YngANDOld_SAME)
        X = ssfSubset(:,PCs_YngANDOld_SAME);
    % Fits a model of ones
    Y = ones(N,1);
    beta = inv(X'*X)*X'*Y;
    SSF_YngANDOld_SAME = X*beta;
end
% Young only
if ~isempty(PCs_YNG)
    X = ssfSubset(:,PCs_YNG);
    % Fits a model of ones
    Y = ones(N,1);
    Y(Oind) = 0;
    beta = inv(X'*X)*X'*Y;
    SSF_YNG = X*beta;
end
% Old only
if ~isempty(PCs_OLD)
    X = ssfSubset(:,PCs_OLD);
    % Fits a model of ones
    Y = ones(N,1);
    Y(Yind) = 0;
    beta = inv(X'*X)*X'*Y;
    SSF_OLD = X*beta;
else
    SSF_OLD = zeros(N,1);
end

corr([SSF_YngANDOld_SAME SSF_YngANDOld_DIFF SSF_YNG SSF_OLD  data.Y data.X])
% Do any of these SSFs relate to performance?

S_YngANDOld_SAME = subfnregstats(data.Y,[data.X SSF_YngANDOld_SAME data.COV]);
S_YngANDOld_SAME.tstat.t
S_YngANDOld_DIFF = subfnregstats(data.Y,[data.X SSF_YngANDOld_DIFF data.COV]);
S_YngANDOld_DIFF.tstat.t
S_YNG = subfnregstats(data.Y,[data.X SSF_YNG data.COV]);
S_YNG.tstat.t
S_ALL = subfnregstats(data.Y,[data.X SSF_YngANDOld_SAME SSF_YngANDOld_DIFF SSF_YNG]);
S_ALL.tstat.t

data_YngANDOld_SAME = data;
data_YngANDOld_SAME.M = SSF_YngANDOld_SAME;
P_YngANDOld_SAME = subfnVoxelWiseProcessBatch(data_YngANDOld_SAME);
P_YngANDOld_SAME{1}.AB1{1}.BCaci

data_YngANDOld_DIFF = data;
data_YngANDOld_DIFF.M = SSF_YngANDOld_DIFF;
P_YngANDOld_DIFF = subfnVoxelWiseProcessBatch(data_YngANDOld_DIFF);
P_YngANDOld_DIFF{1}.AB1{1}.BCaci

data_YNG = data;
data_YNG.M = SSF_YNG;
P_YNG = subfnVoxelWiseProcessBatch(data_YNG);
P_YNG{1}.AB1{1}.BCaci


%% Look for a moderating pattern
data_MOD = data;
data_MOD.M = SSF_YngANDOld_DIFF;
data_MOD.V = ssfSubset(:,[1 2 3]);
data_MOD.ModelNum = '14';
coef = subfnFitRegressModelPCs(data_MOD)
NMed = size(data_MOD.M,2);
NMod = size(data_MOD.V,2);
NCov = size(data_MOD.COV,2);
% coefficient for the constant term in the model
const = coef(1);
% coefficients for the mediator
b = coef(1+1:1+NMed);
% coefficient for the x effect
cP = coef(1+NMed+1);
% coefficient(s) for the moderator
v = coef(1+NMed+1+1:(NMod-1+NMed+1+2));
% coefficients for the interaction term
w = coef(1+NMed+1+NMod+1);
% coefficients for the covariates
b_cov = coef(end-NCov+1:end);

% weighted sum of mediator PCs + X effect + weighted sum of moderator PCs +
% interaction effect
if NCov > 0
    fit = ones(N,1)*const + data_MOD.M*b + data_MOD.X*cP + data_MOD.V*v + ((data_MOD.M*b).*(data_MOD.V*v))*w + data_MOD.COV*b_cov;
else
    fit = ones(N,1)*const + data_MOD.M*b' + data_MOD.X*cP + data_MOD.V*v' + ((data_MOD.M*b').*(data_MOD.V*v'))*w;
end
figure(2)
clf
hold on
plot(data_MOD.Y)
plot(fit,'r')


corr(data.Y,fit).^2
%%





PE_behav_fit_coef = subfnregress(data.Y,[ssfSubset(:,selected_PCs) data.X]);
% create the SSF image
temp = eigenimages_noZeroes(:, selected_PCs) * PE_behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent

PE_behav_fit_composite_PC_image = temp / norm(temp);
clear temp;
%%%%% Obtain SSFs associated with the normalized best linear behavioral-fit image %%%%%
PE_behav_fit_composite_PC_image_ssfs = squeeze(data.M) * PE_behav_fit_composite_PC_image;




selected_PCs = find(combo_matrix(find(AIC == min(AIC)),:))

    X = ssfSubset(:,selected_PCs);
    Y = ones(N,1);
    beta = inv(X'*X)*X'*Y;
    fit = X*beta;
    r = Y - fit;
figure(1)
clf
hold on
plot(Y)
plot(fit,'r')



% The aim is to identify a weighted set of PCs that is:
% used by both age groups to the same degree. 
% used by both age groups but to different degrees
% used by the old but not the young
% used by the young but not the old
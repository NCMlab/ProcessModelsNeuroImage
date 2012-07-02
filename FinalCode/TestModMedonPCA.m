BaseDir = '/share/data/users/js2746_Jason/CogReserveAnalyses'
addpath /share/data/users/ch629/New_Code/
addpath /share/data/users/ch629/OrT_Module/
addpath /share/data/users/js2746_Jason/CogReserveAnalyses/FinalCode/
cd(BaseDir)
% load up the results of a PCA 
load workspace_snapshot_20120523
% load up performance
load sSTM_SSF
load sRET_SSF
load sPRO_SSF
load iSTM_SSF
load iRET_SSF
%load iPRO_SSF
load sRT
load Sex

% matlabpool open 8
%%
N_PC = 12;
NSub = size(vars.ssf,1);
%SSF = sSTM_SSF(:,1:N_PC);
%SSF = sRET_SSF(:,1:N_PC);
%SSF = sPRO_SSF(:,1:N_PC);
%SSF = iSTM_SSF(:,1:N_PC);
SSF = [iSTM_SSF(:,1:6) iRET_SSF(:,1:6)];%  sPRO_SSF(:,1:5)];
Gr = [zeros(75,1); ones(37,1)];

Y = sRT;
M = SSF(:,1:N_PC);
V = SSF(:,1:N_PC);
X = Gr;

STRAT = Gr;

data = {};
data.Y = Y;
data.M = M;
data.X = X;
data.V = [];
data.W = [];
data.Q = [];
data.R = [];
data.Indices=1;
data.COV = Sex;
data.STRAT = STRAT;
data.Thresholds = [0.05 0.01 0.001 0.005];
data.ModelNum = '4';
data.Nboot = 0;

NCov = size(data.COV,2);
% Now this can be done with the other code I have


power_set = boolean_enumeration_f(N_PC);
AIC = inf(2^N_PC-1, 1);
nOfParams = 1;
tic
parfor i=1:2^N_PC-1
    temp = data;
    subset = find(power_set(i, :));
    temp.M = temp.M(:,subset);
    model = [temp.M temp.X temp.COV ones(NSub,1)];
    nOfParams = size(model,2);
    if NSub-nOfParams > 2
        beta = regress(temp.Y,model);
        fit = model*beta;
        AIC(i) = NSub * log(sum((temp.Y-fit).^2) / NSub )  + ...
            2*(nOfParams)*(nOfParams+1)/(NSub-nOfParams-1) + 2*(nOfParams);
    end
end
toc
%%%%% Sort AIC in ascending order %%%%%
[sorted_AIC, pre_sort_index] = sort(AIC);
%%%%% Using results in the original AIC array, create an AIC array for sequential %%%%%
%%%%%   inclusion of predictors starting with the first predictor; sort this AIC  %%%%%
%%%%%   array in a similar fashion as in the original AIC array                   %%%%%
current_row_number = 0;
AIC_sequential_PCs = zeros(N_PC, 1);
for i=(N_PC-1):-1:0
    current_row_number = current_row_number + 2^i;
    AIC_sequential_PCs(N_PC - i) = AIC(current_row_number);
end
[sorted_AIC_sequential_PCs, pre_sort_index_sequential_PCs] = sort(AIC_sequential_PCs);


disp([int2str(2^N_PC-1) ' sets of predictors were tested.'])
disp(' ')

disp('Under arbitrary combination of predictors, the smallest three Akaike Information Criteria values and their respective predictors are as follows:')
for i=1:3
    disp(['     (', int2str(i), ') ', num2str(sorted_AIC(i)), ...
        ' from fitting the following predictors ---> ', ...
        int2str(find(power_set(pre_sort_index(i), :))) ])
end
disp(' ')

disp('Under sequential inclusion of predictors, the smallest three Akaike Information Criteria values and their respective predictors are as follows:')
for i=1:3
    disp(['     (', int2str(i+3), ') ', num2str(sorted_AIC_sequential_PCs(i)), ...
        ' from fitting the following predictors ---> ', ...
        int2str(1:pre_sort_index_sequential_PCs(i)) ])
end
disp(' ')

%
%MEDlist = find(power_set(pre_sort_index(1),:));
MEDlist = 1:pre_sort_index_sequential_PCs(1);
MEDlist = 1:12;
temp = data;
temp.Nboot = 5000;
temp.M = temp.M(:,MEDlist);

[ParameterToBS ] = subfnProcessModelFit(temp,'4',1);

MedSSF = temp.M*ParameterToBS;
temp.M = MedSSF;
Parameters = subfnVoxelWiseProcessBatch(temp);


fprintf(1,'======================================================\n');
fprintf(1,'Indirect effect of X on Y via M\n');
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n','PC#','Effect','Boot SE','BootLLCI','BootUPCI');
for i = 1:length(Parameters{1}.AB)
    fprintf(1,'%10d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',i,Parameters{1}.AB{i}.pointEst,Parameters{1}.AB{i}.se,Parameters{1}.AB{i}.BCaci.alpha05(1),Parameters{1}.AB{i}.BCaci.alpha05(2));
end

fprintf(1,'Effect of X on M\n');
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n','PC#','coeff','SE','t','p');
for i = 1:length(Parameters{1}.AB)
    fprintf(1,'%10d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',i,Parameters{1}.A{i}.beta,Parameters{1}.A{i}.se,Parameters{1}.A{i}.t,Parameters{1}.A{i}.p);
end
   
fprintf(1,'Effect of M on Y\n');
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n','PC#','coeff','SE','t','p');
for i = 1:length(Parameters{1}.AB)
    fprintf(1,'%10d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',i,Parameters{1}.B{i}.beta,Parameters{1}.B{i}.se,Parameters{1}.B{i}.t,Parameters{1}.B{i}.p);
end


%% Now find the moderator
N_PC = 10;
tempMOD = temp;
tempMOD.V = SSF(:,1:N_PC);
N_V = N_PC;

power_set = boolean_enumeration_f(N_PC);
AIC = inf(2^N_PC-1, 1);
nOfParams = 1;
N_M = 1;
N = length(tempMOD.Y);
weights = {};
weights.med = zeros(N_M,1);
weights.mod = zeros(N_V,1);
weights.w = 0;
weights.const = 0;
weights.cP = 0;
W = [zeros(1,size(tempMOD.M,2)+size(tempMOD.V,2)+1+1+1)];
options = optimset('TolX',1e-10,'MaxIter',15000,'MaxFunEvals',15000,'Algorithm','levenberg-marquardt','TolFun',1e-10,'GradObj','on');
tic
parfor i=1:2^N_PC-1
    temp2 = tempMOD;
    subset = find(power_set(i, :));
    temp2.V = temp2.V(:,subset);
    NMod = length(subset);
    nOfParams = 4 + NMod;
    if NSub-nOfParams > 2
        [Wout ,FVAL,EXITFLAG,OUTPUT] = fminsearch('subfnModMedFMIN',W,options,tempMOD);
        b = Wout(1);
        v = Wout(2:1+NMod)';
        w = Wout(1+NMod+1);
        const = Wout(1+NMod+2);
        cP = Wout(1+NMod+3);
        fit = temp2.X*cP + temp2.M*b + temp2.V*v + w*((temp2.M*b).*(temp2.V*v)) + const*ones(N,1);
        AIC(i) = NSub * log(sum((temp2.Y-fit).^2) / NSub )  + ...
            2*(nOfParams)*(nOfParams+1)/(NSub-nOfParams-1) + 2*(nOfParams);
    end
end
toc
%% %%% Sort AIC in ascending order %%%%%
[sorted_AIC, pre_sort_index] = sort(AIC);


%%%%% Using results in the original AIC array, create an AIC array for sequential %%%%%
%%%%%   inclusion of predictors starting with the first predictor; sort this AIC  %%%%%
%%%%%   array in a similar fashion as in the original AIC array                   %%%%%
current_row_number = 0;
AIC_sequential_PCs = zeros(N_PC, 1);
for i=(N_PC-1):-1:0
    current_row_number = current_row_number + 2^i;
    AIC_sequential_PCs(N_PC - i) = AIC(current_row_number);
end
[sorted_AIC_sequential_PCs, pre_sort_index_sequential_PCs] = sort(AIC_sequential_PCs);

disp([int2str(2^N_PC-1) ' sets of predictors were tested.'])
disp(' ')

disp('Under arbitrary combination of predictors, the smallest three Akaike Information Criteria values and their respective predictors are as follows:')
for i=1:3
    disp(['     (', int2str(i), ') ', num2str(sorted_AIC(i)), ...
        ' from fitting the following predictors ---> ', ...
        int2str(find(power_set(pre_sort_index(i), :))) ])
end
disp(' ')

disp('Under sequential inclusion of predictors, the smallest three Akaike Information Criteria values and their respective predictors are as follows:')
for i=1:3
    disp(['     (', int2str(i+3), ') ', num2str(sorted_AIC_sequential_PCs(i)), ...
        ' from fitting the following predictors ---> ', ...
        int2str(1:pre_sort_index_sequential_PCs(i)) ])
end
disp(' ')
MODlist = find(power_set(pre_sort_index(1),:));

%%

tempMOD.V = tempMOD.V(:,MODlist);
[Wout ,FVAL,EXITFLAG,OUTPUT] = fminsearch('subfnModMedFMIN',W,options,tempMOD);



b = Wout(1:1);
v = Wout(1+1:1+N_V)';
w = Wout(1+N_V+1);
const = Wout(1+N_V+2);
cP = W(1+N_V+3);

MOD = v;
tempMOD.V = tempMOD.V*v;

tempMOD.ModelNum = '14'
Parameters = subfnVoxelWiseProcessBatch(tempMOD)




%%
tempr = corr([medSSF modSSF]);
tempr = corr([MedWeights ModWeights]);
% r(i) = tempr(1,2);

% end
% figure
% plot(r)
% create the covariance patterns
tag = 'Gr_sRET_sRT_MED'
selected_PCs = 1:N_PC;
temp = vars.eigenimages_noZeroes(:, selected_PCs) * MedWeights;  %nuisance regressors stay silent
behav_fit_composite_PC_image              = zeros(vars.dim(1:3));
behav_fit_composite_PC_image(vars.meanful_set) = temp / norm(temp);
clear temp;
V = vars.header_info;
V.fname = fullfile(BaseDir,[tag '.nii']);
spm_write_vol(V,behav_fit_composite_PC_image);

tag = 'Gr_sRET_sRT_MOD'
selected_PCs = 1:N_PC;
temp = vars.eigenimages_noZeroes(:, selected_PCs) * ModWeights';  %nuisance regressors stay silent
behav_fit_composite_PC_image              = zeros(vars.dim(1:3));
behav_fit_composite_PC_image(vars.meanful_set) = temp / norm(temp);
clear temp;
V = vars.header_info;
V.fname = fullfile(BaseDir,[tag '.nii']);
spm_write_vol(V,behav_fit_composite_PC_image);

    
fid1 = fopen('medSSF.txt','w');
fid2 = fopen('modSSF.txt','w');
for i = 1:112
    fprintf(fid1,'%0.6f\n',medSSF(i));
    fprintf(fid2,'%0.6f\n',modSSF(i));
end
fclose(fid1)
fclose(fid2)

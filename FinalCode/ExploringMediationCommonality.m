
N = 300;
a = -0.46;
b = 0.03;
cP = -0.42;



Gr = [zeros(N/2,1); ones(N/2,1)];

M = randn(N,1)*0.5 + a.*Gr;

Y = cP.*Gr + randn(N,1)*0.01 + b.*M;
corr([Gr M Y])

%Y = Y - mean(Y);
%M = M - mean(M);

%Y = Y./std(Y);
%M = M./std(M);
%Gr = Gr - mean(Gr);
%Gr = Gr./std(Gr);



data = {};
data.names = {};
data.names.X = 'A';
data.names.M = {'B'};
data.names.Y = 'C';
data.names.V = '';
data.names.W = '';
data.names.Q = '';
data.names.R = '';
data.names.COV = {};
data.Y = Y;
data.M = M;
data.X = Gr;
data.COV = [];
data.Q = [];
data.W = [];
data.R = [];
data.V = [];
data.Indices = 1;
data.STRAT = [];
data.ModelNum = '4';
data.Nboot = 5000;
data.Thresholds = [0.05 0.01 0.005 0.001];
Parameters = subfnVoxelWiseProcessBatch(data);

S = regstats(data.Y,[data.X]);

%
voxel = 1;
fprintf(1,'==========================================================\n');
fprintf(1,'Model Coefficients\n')
fprintf('Independent variable=M\n')
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n',' ','Effect','se','t','p')
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','X',Parameters{voxel}.A{1}.beta,Parameters{voxel}.A{1}.se,Parameters{voxel}.A{1}.t,Parameters{voxel}.A{1}.p)
fprintf(1,'----------------------------------------------------------\n');
fprintf(1,'Model Coefficients\n')
fprintf('Independent variable=Y\n')
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n',' ','Effect','se','t','p')
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','X',Parameters{voxel}.C.beta,Parameters{voxel}.C.se,Parameters{voxel}.C.t,Parameters{voxel}.C.p)
fprintf(1,'----------------------------------------------------------\n');
fprintf(1,'Model Coefficients\n')
fprintf('Independent variable=Y\n')
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n',' ','Effect','se','t','p')
%fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','Constant',Constant,se(1,1),S.tstat.t(1),S.tstat.pval(1))
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','M',Parameters{voxel}.B{1}.beta,Parameters{voxel}.B{1}.se,Parameters{voxel}.B{1}.t,Parameters{voxel}.B{1}.p)
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','X',Parameters{voxel}.CP.beta,Parameters{voxel}.CP.se,Parameters{voxel}.CP.t,Parameters{voxel}.CP.p)
fprintf(1,'----------------------------------------------------------\n');
fprintf(1,'Indirect Effect, alpha = 0.05\n')
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n',' ','Effect','BootSE','BootLower','BootUpper')
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','AB',Parameters{voxel}.AB{1}.pointEst,Parameters{voxel}.AB{1}.se,Parameters{voxel}.AB{1}.BCaci.alpha05(1),Parameters{voxel}.AB{1}.BCaci.alpha05(2));
fprintf(1,'----------------------------------------------------------\n');
fprintf(1,'Model Coefficients\n')
fprintf('Independent variable=Y\n')
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n',' ','Effect','se','t','p')
%fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','Constant',Constant,se(1,1),S.tstat.t(1),S.tstat.pval(1))
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','M',S.beta(2),S.tstat.se(2),S.tstat.t(2),S.tstat.pval(2));

fprintf(1,'----------------------------------------------------------\n');

%%

[C, V, T,r, R2, R2_13, R2_12, R2_23] = subfnCommonality(Y, [Gr M]);
C
V
[error img] = Commonality_2Pred(Y,Gr,M);

%%
N = 100;
voxel = 1;
a = [0:0.01:1];
Nsteps = length(a);
b = 1;
c = 1;
X = [zeros(N/2,1);ones(N/2,1)];
output = zeros(Nsteps,8);
key = {'a','c','b','cP','ab','COM','Ux','Um'};
for i = 1:Nsteps
    M = randn(N,1)*1 + a(i).*X;
    Y = c.*X + randn(N,1)*0.1 + b.*M;
    data = {};
    data.Y = Y;
    data.M = M;
    data.X = X;
    data.COV = [];
    data.Q = [];
    data.W = [];
    data.R = [];
    data.V = [];
    data.Indices = 1;
    data.STRAT = [];
    data.ModelNum = '4';
    data.Nboot = 10;
    data.Thresholds = [0.05 0.01 0.005 0.001];
    
    Parameters = subfnVoxelWiseProcessBatch(data);
    [C, V, T,r, R2, R2_13, R2_12, R2_23] = subfnCommonality(Y, [X M]);
    
    output(i,1) = Parameters{voxel}.A{1}.beta;
    output(i,2) = Parameters{voxel}.C.beta;
    output(i,3) = Parameters{voxel}.B{1}.beta;
    output(i,4) = Parameters{voxel}.CP.beta;
    output(i,5) = Parameters{voxel}.AB{1}.pointEst;
    output(i,6) = sqrt(C);
    output(i,7) = sqrt(V(1,1));
    output(i,8) = sqrt(V(2,2));
end

figure(12)
plot(output(:,1),output(:,5:8))
xlabel('a')
legend('ab','common','cYM','cYX')
figure(13)
plot(output(:,5),output(:,6:8))
xlabel('ab')
legend('common','cYM','cYX')



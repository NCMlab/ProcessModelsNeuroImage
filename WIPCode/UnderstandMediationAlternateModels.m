N = 1000;
A = randn(100,1)*4+0.1;
B = randn(100,1)*2-0.3;% + 0.5.*A;
C = randn(100,1)*3+0.2;% + 0.5.*A;

zA = (A-mean(A))./std(A);
zB = (B-mean(B))./std(B);
zC = (C-mean(C))./std(C);

%[mean(A) std(A) mean(B) std(B)]
%[mean(zA) std(zA) mean(zB) std(zB)]

X = [A B];
zX = [zA zB];


%[A'*A zA'*zA B'*B zB'*zB]
%[X'*X zX'*zX]./N

%inv(A'*A)*A'*B
inv(zA'*zA)*zA'*zC
inv(zX'*zX)*zX'*zC

zX'*zC./N

corr([A B C zA zB zC])
cov([A B C zA zB zC])


%% 
clear
close all
N = 1000;
ABC = randn(N,1);
AB_C = randn(N,1);
BC_A = randn(N,1);
AC_B = randn(N,1);
abc = 0.5;
ab_c = 1;
ac_b = 1;
bc_a = 1;
A = ab_c*AB_C + ac_b*AC_B + abc*ABC + randn(N,1);
B = ab_c*AB_C + bc_a*BC_A + abc*ABC + randn(N,1);
C = ac_b*AC_B + bc_a*BC_A + abc*ABC + randn(N,1);
A = (A - mean(A))./std(A);
B = (B - mean(B))./std(B);
C = (C - mean(C))./std(C);

corr([A B C])
S1 = regstats(C,[A B]);
S1c= regstats(C,[A]);
%S1a = regstats(B,A);
eff1 = [S1.beta(2) S1c.beta(2) S1c.beta(2) - S1.beta(2)];

S2 = regstats(B,[A C]);
S2c = regstats(B,[A]);
eff2 = [S2.beta(2) S2c.beta(2) S2c.beta(2) - S2.beta(2)];

S3 = regstats(C,[B A]);
S3c = regstats(C,[B]);
eff3 = [S3.beta(2) S3c.beta(2) S3c.beta(2) - S3.beta(2)];



names = {'A' 'B' 'C'};
[error img] = Commonality_2Pred(C,A,B,names);

%%
% run mediation analyses 
%
data = {};
data.Xname = 'A';
data.Yname = 'C';
data.Mname = 'B';
data.Vname = 'V';
data.Wname = 'W';

data.STRAT = [];
data.COV = [];%randn(N,2);
data.V = [];
data.W = [];
data.Q = [];
data.R = [];
data.ModelNum = '4';
data.Thresholds = [0.05 0.01 0.005];
data.Indices = 1;
data.Nboot = 2000;

% Calculate the full stats of the model
%[ParameterToBS Parameters] = subfnProcessModelFit(data,data.ModelNum,PointEst);

%tic
data.Xname = 'A';
data.Yname = 'C';
data.Mname = 'B';
data.X = A;
data.Y = C;
data.M = B;
Parameters_AC_B = subfnVoxelWiseProcessBatch(data);

data.Xname = 'A';
data.Yname = 'B';
data.Mname = 'C';
data.X = A;
data.Y = B;
data.M = C;
Parameters_AB_C = subfnVoxelWiseProcessBatch(data);

data.Xname = 'B';
data.Yname = 'C';
data.Mname = 'A';
data.X = B;
data.Y = C;
data.M = A;
Parameters_BC_A = subfnVoxelWiseProcessBatch(data);


fprintf(1,'=======================================\n');
fprintf(1,'%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\n','Model','TotalEff','DirectEff','DiffEff',...
    'IndirectEff','lowCI','hiCI');
fprintf(1,'%8s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n','A-C.B',eff1(2),eff1(1),eff1(3),...
    Parameters_AC_B{1}.AB1{1}.pointEst,...
    Parameters_AC_B{1}.AB1{1}.BCaci.alpha05(1),Parameters_AC_B{1}.AB1{1}.BCaci.alpha05(2));
fprintf(1,'%8s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n','A-B.C',eff2(2),eff2(1),eff2(3),...
    Parameters_AB_C{1}.AB1{1}.pointEst,...
    Parameters_AB_C{1}.AB1{1}.BCaci.alpha05(1),Parameters_AB_C{1}.AB1{1}.BCaci.alpha05(2));
fprintf(1,'%8s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n','B-C.A',eff3(2),eff3(1),eff3(3),...
    Parameters_BC_A{1}.AB1{1}.pointEst,...
    Parameters_BC_A{1}.AB1{1}.BCaci.alpha05(1),Parameters_BC_A{1}.AB1{1}.BCaci.alpha05(2));


% % 
clear

% Test new mediation code
N = 112;
Gr = round(rand(N,1));
Nmed = 1;
V = randn(N,1); 
W = randn(N,1); 

X = randn(N,Nmed) + 10; 
M1 = 0.2*X + randn(N,1)*0.5;% + i;
M2 = 0.2*X + randn(N,1)*0.5;% + i;
Y = 0.3*M1(:,1) + 0*M2(:,1) + 0.5.*X + randn(N,1)*0.5;% + 0.1.*X.*M1(:,1);

% X = randn(N,1);
% M = randn(N,1);
% Y = randn(N,1);

%corr([X M V W Y])
%corr([X M Y])
%regress(Y,[X ones(N,1)])
%regress(M,[X ones(N,1)])
%
% clear
% load ../PracticeData/ModMeddata

data = {};
data.names = {};
data.names.X = 'A';
data.names.M = {'B1' 'B2'};
data.names.Y = 'C';
data.names.V = '';
data.names.W = '';
data.names.Q = '';
data.names.R = '';
data.names.COV = {};




data.STRAT = [];
data.COV = [];%randn(N,2);
data.W = [];
data.Q = [];
data.R = [];
data.V = [];

data.ModelNum = '4';
data.Thresholds = [0.05];
data.Indices = 1;
data.Nboot = 2000;

% Model of interest
data.X = X;
data.Y = Y;
data.M = [M1];

% 
%  data.X = randn(N,1);
%  data.M = randn(N,1) + 0.4.*data.X;
%  data.Y = randn(N,1) + 0.1.*data.M + 0.5.*data.M.*data.X;
 

%data.V = V;
% Add the bootstrap resamples to the input data
NSub = length(data.X);
data2 = data;
%data.Resamples = subfuncBootStrapResamples(data.Nboot,NSub,data.STRAT);
tic
Parameters2 = subfnVoxelWiseProcessBatch(data2);
toc
subfnPrintResults(Parameters2{1})
%% BOOT

Ape = regress(data2.M,[ones(N,1) data2.X]);
Bpe = regress(data2.Y,[ones(N,1) data2.M data2.X]);
ABpe = Ape(2)*Bpe(2);

Nboot = 5000;
BOOTpe = zeros(Nboot,1);
for i = 1:Nboot
    tempData = data2;
    Samp = WIPBootStrap([],N);
    tempData.X = data2.X(Samp);
    tempData.M = data2.M(Samp);
    tempData.Y = data2.Y(Samp);
    A = regress(tempData.M,[ones(N,1) tempData.X]);
    B = regress(tempData.Y,[ones(N,1) tempData.M tempData.X]);
    C = regress(tempData.Y,[ones(N,1) tempData.X]);
    %[pointEst] = subfnProcessModelFit(tempData,0);
    BOOTpe(i) = A(2)*B(2);
end
figure(1)
clf
hist(BOOTpe,40)
%pe = Parameters2{1}.AB1{1}.pointEst;
line([ABpe ABpe],[0 50])

corr([data2.X data2.M data2.Y])
% PERM
Nperm = 5000;
PERMpe = zeros(Nperm,1);
for i = 1:Nperm
    tempData = data2;
    tempData.Nboot = 0;
    tempData.Y = data2.Y(randperm(N));
    A = regress(tempData.M,[ones(N,1) tempData.X]);
    B = regress(tempData.Y,[ones(N,1) tempData.M tempData.X]);
    C = regress(tempData.Y,[ones(N,1) tempData.X]);
    %[pointEst] = subfnProcessModelFit(tempData,0);
    PERMpe(i) = C(2) - B(3);
    PERMpe2(i) = A(2)*B(2);%pointEst.values;
end
figure(2)
hist(PERMpe,40)
%pe = Parameters2{1}.AB1{1}.pointEst;
line([ABpe ABpe],[0 50])
sum(abs(PERMpe)>abs(ABpe))/Nperm

%%
data = [];
data.data{1} = X;
data.data{2} = M;
data.data{3} = Y;
data.Direct = zeros(3,3);
data.Direct(1,[2 3]) = 1;
data.Direct(2,3) = 1;
data.Inter = zeros(3,3);
data.Inter([1 3],3) = 1;
data.Paths = zeros(3,3);
data.Paths(1,2) = 1;
data.Paths(2,3) = 2;
data.Thresholds = data2.Thresholds;
data.Nboot = data2.Nboot;
data.STRAT = data2.STRAT;
data.Indices = 1;
data.NSub = N;
data.Nvar =3 ;
% Create the bootstrap resamples
Samp = int16(zeros(data.NSub,data.Nboot));
for i = 1:data.Nboot;
    Samp(:,i) = WIPBootStrap(data.STRAT,data.NSub);
end
data.BootSamp = Samp;

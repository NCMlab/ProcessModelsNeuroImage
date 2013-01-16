
% functions
% bootSE
% subfnBootStrp
% subfnFindConfidenceIntervals
% subfnProcessModelFit


clear
% Test new mediation code
N = 112;
AonB = 0.35;
BonC = 0.35;
AonC = 0.5;
sC = 0.5;
sB = 0.5;

A = round(rand(N,1));
B = A.*AonB + randn(N,1).*sB;
C = A.*AonC + B.*BonC + randn(N,1).*sC;
corr([A B C])

%Y = 0.25*M(:,1) + randn(N,1)*0.15 + 0.5.*X.*M(:,1);

V = randn(N,1); 
W = randn(N,1); 


% % functions
% % bootSE
% % subfnBootStrp
% % subfnFindConfidenceIntervals
% % subfnProcessModelFit
% 
% 
clear

% Test new mediation code
N = 112;
Gr = round(rand(N,1));
Nmed = 3;
V = randn(N,1); 
W = randn(N,1); 

M = randn(N,Nmed) + 10; 
X = zeros(N,1);
for i = 1:Nmed
    X = X + 0.25*M(:,i) + randn(N,1)*0.15 + i;
end

Y = 0.25*M(:,1) + randn(N,1)*0.15 + 0.5.*X.*M(:,1);

%corr([X M V W Y])
%corr([X M Y])
%regress(Y,[X ones(N,1)])
%regress(M,[X ones(N,1)])
%
% clear
% load ../PracticeData/ModMeddata

data = {};
data.names = {};
data.names.X = 'X';
data.names.M = {'M'};
data.names.Y = 'Y';
data.names.V = 'V';
data.STRAT = [];
data.COV = [];%randn(N,2);
data.V = V;
data.W = V;
data.Q = [];
data.R = [];
data.ModelNum = '74';
data.Thresholds = [0.95];
data.Indices = 1;
data.Nboot = 2000;

% Model of interest

data.X = A;
data.Y = C;
data.M = B;
Parameters = subfnVoxelWiseProcessBatch(data);
subfnPrintResults(Parameters{1})




% % alternate model 1
data.Xname = 'A';
data.Yname = 'B';
data.Mname = 'C';
data.X = A;
data.Y = B;
data.M = C;
Parameters = subfnVoxelWiseProcessBatch(data);
subfnPrintResults(Parameters{1})
% % alternate model 2
data.Xname = 'B';
data.Yname = 'C';
data.Mname = 'A';
data.X = B;
data.Y = C;
data.M = A;
Parameters = subfnVoxelWiseProcessBatch(data);
subfnPrintResults(Parameters{1})
% % alternate model 3
data.Xname = 'C';
data.Yname = 'A';
data.Mname = 'B';
data.X = C;
data.Y = A;
data.M = B;
Parameters = subfnVoxelWiseProcessBatch(data);
subfnPrintResults(Parameters{1})
% % alternate model 4
data.Xname = 'C';
data.Yname = 'B';
data.Mname = 'A';
data.X = C;
data.Y = B;
data.M = A;
Parameters = subfnVoxelWiseProcessBatch(data);
subfnPrintResults(Parameters{1})
% % alternate model 5
data.Xname = 'B';
data.Yname = 'A';
data.Mname = 'C';
data.X = B;
data.Y = A;
data.M = C;

data.Nboot = 100;



% Calculate the full stats of the model
%[ParameterToBS Parameters] = subfnProcessModelFit(data,data.ModelNum,PointEst);

%tic
>>>>>>> 6b66aa81ba5272e82df9b685933e24244a7e1f4a
Parameters = subfnVoxelWiseProcessBatch(data);
subfnPrintResults(Parameters{1})

%%
% plot probe
% Nprobe = length(Parameters{1}.CondAB{1});
% probeValues = zeros(Nprobe-1,1);
% EffectValues  = zeros(Nprobe-1,1);
% se  = zeros(Nprobe-1,1);
% lowerCI = zeros(Nprobe-1,1);
% upperCI = zeros(Nprobe-1,1);
% for i = 2:Nprobe
%     probeValues(i-1) = Parameters{1}.CondAB{1}{i}.probeValue;
%     EffectValues(i-1) = Parameters{1}.CondAB{1}{i}.pointEst;
%     se(i-1) = Parameters{1}.CondAB{1}{i}.bootSE;
%     temp = Parameters{1}.CondAB{1}{i};
%     lowerCI(i-1) = temp.BCaci.alpha05(1);
%     upperCI(i-1) = temp.BCaci.alpha05(2);
% end
% 
% for i = 1:Nprobe - 1 
%     fprintf(1,'%10.4f\t%10.4f\t%10.5f\t%10.4f\t%10.4f\n',probeValues(i),EffectValues(i),se(i),lowerCI(i),upperCI(i));
% end


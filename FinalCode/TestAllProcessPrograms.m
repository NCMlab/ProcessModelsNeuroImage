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
data.Xname = 'X';
data.Yname = 'Y';
data.Mname = 'M';
data.Vname = 'V';
data.Wname = 'W';
data.X = X;
data.Y = Y;
data.M = M;
data.STRAT = [];
data.COV = [];%randn(N,2);
data.V = V;
data.W = V;
data.Q = [];
data.R = [];
data.ModelNum = '14';
data.Thresholds = [0.95];
data.Indices = 1;
data.Nboot = 1000;

% Calculate the full stats of the model
%[ParameterToBS Parameters] = subfnProcessModelFit(data,data.ModelNum,PointEst);

%tic
Parameters = subfnVoxelWiseProcessBatch(data);
%toc

%Parameters{1}
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


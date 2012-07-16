% functions
% bootSE
% subfnBootStrp
% subfnFindConfidenceIntervals
% subfnProcessModelFit


clear
% Test new mediation code
N = 112;
Gr = round(rand(N,1));
Nmed = 3;
M = randn(N,Nmed) + 10; 
X = zeros(N,1);
for i = 1:Nmed
    X = X + 0.25*M(:,i) + randn(N,1)*0.15 + i;
end
Y = 0.25*M(:,1) + randn(N,1)*0.15 + 0.5.*X.*M(:,1);

V = randn(N,1); 
W = randn(N,1); 
%corr([X M V W Y])
%corr([X M Y])
%regress(Y,[X ones(N,1)])
%regress(M,[X ones(N,1)])
%%
data = {};
data.Xname = 'A';
data.Yname = 'B';
data.Mname = 'C';
data.Vname = 'V';
data.Wname = 'W';
data.X = X;
data.Y = M;
data.M = Y;
data.STRAT = Gr;
data.COV = [];
data.V = V;
data.W = V;
data.Q = [];
data.R = [];
data.ModelNum = '4';
data.Thresholds = [0.05 0.01 0.005];
data.Indices = 1;
data.Nboot = 5000;

% Calculate the full stats of the model
%[ParameterToBS Parameters] = subfnProcessModelFit(data,data.ModelNum,PointEst);

%tic
Parameters = subfnVoxelWiseProcessBatch(data);
%toc
%Parameters{1}

%
if ~isfield(data,'Xname')
    data.Xname = 'X';
end
if ~isfield(data,'Yname')
    data.Yname = 'Y';
end
if ~isfield(data,'Mname')
    data.Mname = 'M';
end
if ~isfield(data,'Vname')
    data.Vname = 'V';
end
%%
fprintf(1,'======================================================\n');
fprintf(1,'Model = %s\n',data.ModelNum);
fprintf(1,'\tY = %s\n',data.Yname);
fprintf(1,'\tX = %s\n',data.Xname);
fprintf(1,'\tM = %s\n',data.Mname);

fprintf(1,'Sample size = %d\n\n',length(data.X));

fprintf(1,'Indirect effect of %s on %s via %s (a*b pathway)\n',data.Xname,data.Yname,data.Mname);
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n','PC#','Effect','Boot SE','BootLLCI','BootUPCI');
fprintf(1,'%10d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',i,Parameters{1}.AB1{1}.pointEst,Parameters{1}.AB1{1}.bootSE,Parameters{1}.AB1{1}.BCaci.alpha05(1),Parameters{1}.AB1{1}.BCaci.alpha05(2));

fprintf(1,'Total effect of %s on %s (c pathway)\n',data.Xname,data.Yname);
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n','factor','beta','SE','t','p');
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','const',Parameters{1}.Cconst.beta,Parameters{1}.Cconst.se,Parameters{1}.Cconst.t,Parameters{1}.Cconst.p);
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',data.Xname,Parameters{1}.C.beta,Parameters{1}.C.se,Parameters{1}.C.t,Parameters{1}.C.p);

fprintf(1,'\nEffect of %s on %s (a pathway)\n',data.Xname,data.Mname);
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n','factor','coeff','SE','t','p');
for i = 1:length(Parameters{1}.AB1)
    fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','Const',Parameters{1}.Aconst{i}.beta,Parameters{1}.Aconst{i}.se,Parameters{1}.Aconst{i}.t,Parameters{1}.Aconst{i}.p);
    fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',data.Xname,Parameters{1}.A{i}.beta,Parameters{1}.A{i}.se,Parameters{1}.A{i}.t,Parameters{1}.A{i}.p);
end
   
fprintf(1,'\nEffect of %s on %s (b and c-prime effects)\n',data.Mname,data.Yname);
fprintf(1,'%10s\t%10s\t%10s\t%10s\t%10s\n','factor','coeff','SE','t','p');
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n','Const',Parameters{1}.Bconst.beta,Parameters{1}.Bconst.se,Parameters{1}.Bconst.t,Parameters{1}.Bconst.p);
fprintf(1,'%10s%d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',data.Mname,i,Parameters{1}.B{i}.beta,Parameters{1}.B{i}.se,Parameters{1}.B{i}.t,Parameters{1}.B{i}.p);
fprintf(1,'%10s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',data.Xname,Parameters{1}.CP.beta,Parameters{1}.CP.se,Parameters{1}.CP.t,Parameters{1}.CP.p);

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


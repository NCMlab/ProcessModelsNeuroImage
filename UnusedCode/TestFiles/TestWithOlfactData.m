    
File = 'OlfactData.txt'

[A B C D] = textread(File,'%s%s%s%s','delimiter',',');
X = str2num(char(A)); % group
Y = str2num(char(C)); % perf
M = str2num(char(B)); % thresh
V = str2num(char(D)); % piriform
data = {};
data.X = X;
data.Y = Y;
data.M = M;
data.V = V;
data.W = [];
data.Q = [];
data.R = [];
data.COV = [];
data.STRAT = [];

ModelNum = '14';
Thresholds = [0.05 0.01 0.005 0.001];
Nboot = 1000;
data.ModelNum = ModelNum;
data.Thresholds = Thresholds;
data.Nboot = Nboot;
Parameters = subfnVoxelWiseProcessBatch(data)



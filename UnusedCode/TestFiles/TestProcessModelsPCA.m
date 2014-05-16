% function TestProcessModelsPCA
clear
cd('C:\Users\Lloyd Hastings\Desktop\SteffenerColumbia\AcceptedPapers\NeuralSignaturesOfCognition\Data')
Data = textread('101712.csv','%s','delimiter',',');

NSub = 58;
Data = reshape(Data,length(Data)/(NSub+1),(NSub+1))';

AgeGroup = str2num(char(Data(2:end,3)));
education = str2num(char(Data(2:end,4)));
NARTerr = str2num(char(Data(2:end,9)));
MEM = str2num(char(Data(2:end,39)));
CBF = reshape(str2num(char(Data(2:end,108:119))),NSub,12);

% Pick the best combo of CBF 
M = CBF(:,[2 7 8 12]);
S = subfnregstats(MEM,[M AgeGroup])

CBFmem = M*S.beta(2:5);


data = {};
data.Xname = 'Age';
data.Yname = 'Mem';
data.Mname = 'CBF';
data.Vname = 'NARTerr';
data.Wname = 'NARTerr';
data.X = AgeGroup;
data.Y = MEM;
data.M = CBF;
data.STRAT = [];
data.COV = [];%randn(N,2);
data.V = NARTerr;
data.W = NARTerr;
data.Q = [];
data.R = [];
data.ModelNum = '58';
data.Thresholds = [0.05 0.01];
data.Indices = 1;
data.Nboot = 100;



clear
N = 100;
Nmed = 4;
data = {};
data.X = round(rand(N,1));
AgeEffectsOnM = rand(Nmed,1);

data.V = randn(N,1);

data.M = zeros(N,Nmed);
for i = 1:Nmed
    data.M(:,i) = data.X*AgeEffectsOnM(i) + 0.1.*randn(N,1);
end
MeffectOnY = [0.8 0 0 3]';
cPIn = 0.2;
Vin = 0;
InterIn = 0;
ConstIn = 10;
Win = [MeffectOnY; Vin; InterIn; ConstIn; cPIn];
data.Y = data.M*MeffectOnY + data.X*cPIn + data.V*Vin + (InterIn*(data.M*MeffectOnY.*data.V)) + ones(N,1)*ConstIn;
data.Y = data.Y + 0.2.*randn(N,1);
data.STRAT = [];
data.COV = [];
data.Thresholds = [0.05];
%corr([data.X data.M data.V data.Y])

% get data sizes
N = length(data.Y);
NMed = size(data.M,2);
NMod = size(data.V,2);
% initialized parameter vector
W = zeros(NMed + NMod + 3,1);
options = optimset('TolX',1e-10);
options = optimset(options,'MaxIter',15000);
options = optimset(options,'MaxFunEvals',15000);
options = optimset(options,'Algorithm','levenberg-marquardt');
options = optimset(options,'TolFun',1e-10);
options = optimset(options,'GradObj','on');


[Wout ,FVAL,EXITFLAG,OUTPUT] = fminsearch('subfnRegressPCs',W,options,data);

[Win Wout]
% coefficients for the mediator
b = Wout(1:NMed);
% coefficient(s) for the moderator
v = Wout(NMed+1:NMed+NMod);
% coefficients for the interaction term
w = Wout(NMed+NMod+1);
% coefficient for the constant term in the model
const = Wout(NMed+NMod+2);
% coefficient for the x effect
cP = Wout(NMed+NMod+3);

fit = data.M*b + data.X*cP + data.V*v + ((data.M*b).*(data.V*v))*w + ones(N,1)*const;
resid = data.Y - fit;
NParam = length(Win);
NSub = N;

AIC = NSub * log(resid'*resid / NSub )  + ...
    2*(NParam)*(NParam+1)/(NSub-NParam-1) + 2*(NParam);


 Nboot = 500;
 Wboot = zeros(length(W),Nboot);
 for i = 1:Nboot
     Wboot(:,i) = fminsearch('subfnRegressPCs',W,options,subfnReturnBootStrapData(data));
 end
 %%
 TestParam = 1;
alpha = 0.05;
bstat = Wboot(TestParam,:)';
Sbstat = sort(bstat);
i=1;j=1;k=1;
PERci = [Sbstat(ceil(data.Thresholds(i)/2*Nboot),j,k) Sbstat(ceil((1-data.Thresholds(i)/2)*Nboot),j,k)];
[Win(TestParam) Wout(TestParam) PERci]

figure(10)
hist(Wboot(TestParam,:))
%%

figure(1)
clf
hold on
plot(data.Y);
plot(fit,'r')
corrcoef([data.Y fit]).^2


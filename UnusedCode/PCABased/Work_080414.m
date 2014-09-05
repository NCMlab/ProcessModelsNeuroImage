%%
clear
N = 150;
data = {};
data.ModelNum = '14';
data.X = round(rand(N,1));
Npcs = 4;
data.M = zeros(N,Npcs);
for i = 1:Npcs
    data.M(:,i) = data.X.*0.5 + randn(N,1);
end
data.V = randn(N,Npcs);
b_Vpc = rand(Npcs,1);
b_Mpc = rand(Npcs,1);

data.Y = data.M*b_Mpc.*0.5 + randn(N,1) + data.V*b_Vpc;
data.COV=[];

[beta, FVAL, EXITFLAG, fit] = subfnCallRegressPCs(data,beta+randn(size(beta)).*0.1);
FVAL

beta
%%
figure(10)
clf
hold on
plot(data.Y)
plot(fit.Overall,'r')
%%
S1 = regstats(data.Y,[data.M data.X]);
S1.beta


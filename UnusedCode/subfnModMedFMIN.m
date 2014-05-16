function F = subfnModMedFMIN(W, data)


NMed = size(data.M,2);
NMod = size(data.V,2);
b = W(1:NMed);
v = W(NMed+1:NMed+NMod)';
w = W(NMed+NMod+1);
const = W(NMed+NMod+2);
cP = W(NMed+NMod+3);

 
Y = data.Y;
M = data.M;
V = data.V;
N = length(data.Y);
X = data.X;

est = X*cP + M*b + V*v + w*((M*b).*(V*v)) + const*ones(N,1);
%est = X*cP + M*b + const*ones(N,1);

err = (Y - est);
F = err'*err;


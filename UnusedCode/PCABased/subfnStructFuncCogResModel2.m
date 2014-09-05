function F = subfnStructFuncCogResModel2(coef,data)
% P <- A*b_A2 + S*b_S2 + CR*b_CR2 + (CR*S)*b_CRxS2 + Sum(F_pc*b_Fpc)

N = size(data,1);
N_Spc = data(1,end);
N_Fpc = data(2,end);

A = data(:,2);
S = data(:,3);
CR = data(:,4);
estF =  data(:,5:5+N_Fpc-1)*coef(6:6+N_Fpc-1)';
P = data(:,5+N_Fpc);

% P = data(:,end-1);

fit = ones(N,1)*coef(1) + ...
    A*coef(2) + ... % Age
    S*coef(3) + ...% Structure
    CR*coef(4)+ ... % CogRes
    (S.*CR)*coef(5) + ...% Function
    estF;
    %COV*coef(2+S_Npc+F_Npc+3:2+S_Npc+F_Npc+2+COV_N);
    
err = (P - fit);
F = err'*err;
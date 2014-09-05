function fit = subfnStructFuncCogResModelFIT(coef,data)
N = size(data,1);
S_Npc = coef(1);
F_Npc = coef(2);
COV_N = coef(3);

S = data(:,6:6+S_Npc-1)*coef(6:6+S_Npc-1);
demeanS = S - mean(S);

CR = data(:,5+S_Npc+F_Npc+1);
demeanCR = CR - mean(CR);

SxCR = demeanCR.*demeanS;

P = data(:,5+S_Npc+F_Npc+2);
COV = data(:,7+S_Npc + F_Npc + 1: 7+S_Npc + F_Npc + COV_N);


fit = ones(N,1)*coef(3) + ...
    data(:,4)*coef(4) + ... % Age
    S + ...% Structure
    CR*coef(4+S_Npc+F_Npc+1) + ... % CogRes
    SxCR*coef(4+S_Npc+F_Npc+2) + ... % Interaction
    COV*coef(7+S_Npc + F_Npc + 1: 7+S_Npc + F_Npc + COV_N);
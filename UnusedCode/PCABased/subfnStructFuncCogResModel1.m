function F = subfnStructFuncCogResModel1(coef,data)
% P <- A*b_A1 + Sum(S_pc*b_Spc) + CR*b_CR1 + (CR*Sum(S_pc*b_Spc))*b_CRxS1

N = size(data,1);
S_Npc = data(1,end);


A = data(:,2);

S = data(:,3:3+S_Npc-1);
estS = S*coef(3:3+S_Npc-1)';

% demeanS = S - mean(S);
% F = data(:,3+S_Npc:3+S_Npc+F_Npc-1);
% 
CR = data(:,2+S_Npc+1);
%demeanCR = CR - mean(CR);
% 
 %SxCR = demeanCR.*demeanS;
% 
P = data(:,2+S_Npc+2);

% P = data(:,end-1);

fit = ones(N,1)*coef(1) + ...
    A*coef(2) + ... % Age
    estS + ...% Structure
    CR*coef(2+S_Npc+1) + ...% + ... % CogRes
    (estS.*CR)*coef(2+S_Npc + 2);
    %SxCR*coef(2+S_Npc+F_Npc+2) + ... % Interaction
    %COV*coef(2+S_Npc+F_Npc+3:2+S_Npc+F_Npc+2+COV_N);
    
err = (P - fit);
F = err'*err;
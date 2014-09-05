function F = subfnStructFuncCogResModel3(coef,data)
% F <- (CR*S)*b_CRxS3 + CR*b_CR3 + S*b_S3 + A*b_A3

N = size(data,1);
S_Npc = data(1,end);
F_Npc = data(2,end);

A = data(:,2);

S = data(:,3:3+S_Npc-1);

estS = S*coef(3:3+S_Npc-1);
% demeanS = S - mean(S);
% F = data(:,3+S_Npc:3+S_Npc+F_Npc-1);
% 
CR = data(:,2+S_Npc+F_Npc+1);
%demeanCR = CR - mean(CR);
% 
 %SxCR = demeanCR.*demeanS;
% 
P = data(:,2+S_Npc+F_Npc+2);
COV = data(:,2+S_Npc+F_Npc+3:2+S_Npc+F_Npc+2+COV_N);
% P = data(:,end-1);

fit = ones(N,1)*coef(1) + ...
    A*coef(2) + ... % Age
    estS + ...% Structure
    CR*coef(2+S_Npc+F_Npc+1);% + ... % CogRes
 %*coef(3+S_Npc:3+S_Npc+F_Npc-1) + ...% Function
    %SxCR*coef(2+S_Npc+F_Npc+2) + ... % Interaction
    %COV*coef(2+S_Npc+F_Npc+3:2+S_Npc+F_Npc+2+COV_N);
    
err = (P - fit);
F = err'*err;
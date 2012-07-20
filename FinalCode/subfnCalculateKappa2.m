
function k2 = subfnCalculateKappa2(X, M, Y, a, b)
S = cov([X M Y]);
%beta = subfnregress(C,[A B])
%b = beta(3);
%beta = subfnregress(B,[A])
%a = beta(2);
perm_a_1 = sqrt(S(2,3))*sqrt(S(1,3));
perm_a_2 = sqrt(S(2,2)*S(3,3) - S(3,2));
perm_a_3 = sqrt(S(1,1)*S(3,3) - S(3,1));
perm_a_4 = S(1,1)*S(3,3);
perm_a = [(perm_a_1 - perm_a_2*perm_a_3)/perm_a_4 (perm_a_1 + perm_a_2*perm_a_3)/perm_a_4];
if sign(a) > 0
    max_a = perm_a(2);
else
    max_a = perm_a(1);
end
perm_b_1 = sqrt(S(1,1)*S(3,3) - S(3,1));
perm_b_2 = sqrt(S(1,1)*S(2,2) - S(2,1));
perm_b = [-perm_b_1/perm_b_2 perm_b_1/perm_b_2];
if sign(b) > 0
    max_b = perm_b(2);
else
    max_b = perm_b(1);
end

perm_ab = max_a*max_b;
k2 = (a*b)/perm_ab;
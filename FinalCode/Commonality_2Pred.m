
function [error img] = Commonality_2Pred(Y,X1,X2,names)


[B] = regress(Y,[X1 X2 ones(size(X1))]);
pred = [X1 X2 ones(size(X1))]*B;
R2_YX1X2 = corr(Y, pred)^2;

R2_X1X2 = corr(X1, X2).^2;


R2_YX1 = corr(Y,X1)^2;
R2_YX2 = corr(Y,X2)^2;

UX1 = R2_YX1X2 - R2_YX2;
UX2 = R2_YX1X2 - R2_YX1;
CYX1X2 = R2_YX1 + R2_YX2 - R2_YX1X2;

JX1X2_NoY = R2_X1X2 - CYX1X2;

UY = 1 - CYX1X2 - UX1 - UX2;

%               |A|
%               |A and B|
%               |B|
%               |B and C|
%               |C|
%               |C and A|
%               |A and B and C|
%        

[error img] = vennX([UY UX1 [1-CYX1X2-JX1X2_NoY-UX1] [JX1X2_NoY] [1-CYX1X2-JX1X2_NoY-UX2] UX2 CYX1X2],0.005,names);


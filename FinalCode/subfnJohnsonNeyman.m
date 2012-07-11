function JNvalue = subfnJohnsonNeyman(maineffW,maineffCov, intW, intCov,meIntCov, tcrit)
% assuming one mediator
% if the returned values are REAL then there is a range of values where the
% interaction is significant.
% If the returned values are COMPLEX then there is no range of
% significance.
% 
JNvalue = -99;
% 
b1 = maineffW; % mediator weight
%b2 = Stats.beta(3);
b3 = intW; % interaction weight
sb1b3 = meIntCov;

s2b1 = maineffCov;
s2b3 = intCov;
A = -2*((tcrit^2)*sb1b3 - b1*b3);
B = (2*(tcrit^2)*sb1b3 - 2*b1*b3)^2;
C = -4*((tcrit^2)*s2b3 - b3^2)*((tcrit^2)*s2b1 - b1^2);
D = 2*((tcrit^2)*s2b3 - b3^2);

tempJN = [(A + sqrt(B + C))/D (A - sqrt(B+C))/D];
% check to see if there is a signifcant range, otherwise set to -99
if isreal(tempJN)
    JNvalue = tempJN;
end
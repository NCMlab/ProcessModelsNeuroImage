function JNvalue = JohnsonNeyman(Stats,probeInt,tcrit)
b1 = Stats.beta(2); % Main effect of regressor 1
%b2 = Stats.beta(3); % Main effect of regressor 2
b3 = Stats.beta(4); % Interaction between regressors 1 and 2
sb1b3 = Stats.covb(2,4);


s2b1 = Stats.covb(2,2);
s2b3 = Stats.covb(4,4);
A = -2*((tcrit^2)*sb1b3 - b1*b3);
B = (2*(tcrit^2)*sb1b3 - 2*b1*b3)^2;
C = -4*((tcrit^2)*s2b3 - b3^2)*((tcrit^2)*s2b1 - b1^2);
D = 2*((tcrit^2)*s2b3 - b3^2);

JNvalue = [(A + sqrt(B + C))/D (A - sqrt(B+C))/D];
% find the J-N values within the range of M values
outJN = [];
for i = 1:2
    if (JNvalue(i) < probeInt.probeM(end)) & (JNvalue(i) > probeInt.probeM(1))
        outJN = [outJN; JNvalue(i)];
    end
end
JNvalue = outJN;



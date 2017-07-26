function [Z, p, BCaCI] = CalculateBCaCIsimple(PointEstimates,BootStraps, JackKnifes, alpha)
[m, n] = size(PointEstimates);
Z = zeros(m,n);
p = zeros(m,n);
if (m > 1) || (n > 1)
    for i = 1:m
        for j = 1:n
            if PointEstimates(i,j) ~= 0
                [Alpha1, Alpha2, Z(i,j), p(i,j), BCaCI] = CalculateBCaLimitsOneValue(squeeze(JackKnifes(i,j,:)),PointEstimates(i,j), squeeze(BootStraps(i,j,:)),alpha);
            end
        end
    end
else
    [Alpha1, Alpha2, Z p, BCaCI] = CalculateBCaLimitsOneValue(JackKnifes,PointEstimates, BootStraps,alpha);
end
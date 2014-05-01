function BCaCI = WIPsubfnCalculateBCaCI(data,Alpha1,Alpha2,PointEstimate)
% calculate the bias corrected accelerated confidence intervals from the
% pre-calculated alpha limites.

[m n] = size(PointEstimate);
BCaCI = zeros([size(Alpha1) 2]);
Nboot = size(data,3);
for i = 1:m
    for j = 1:n
        if PointEstimate(i,j) ~= 0
            SortBoot = squeeze(sort(data(i,j,:)));
            BCaCI(i,j,:) = [SortBoot(ceil((Alpha1(i,j)*Nboot))) SortBoot(ceil((Alpha2(i,j)*Nboot)))];
        end
    end
end

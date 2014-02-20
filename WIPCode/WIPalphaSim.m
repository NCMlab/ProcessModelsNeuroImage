% THe idea is to create random fields, then smooth them.
m = 100;
n = 100;

a = randn(m,n);

b = randn(m,n);

c = randn(m,n);


alpha = 0.005;
% multiply them
ab = a.*b;
abc = a.*b.*c;
% threshold them
Sabc = sort(reshape(abc,m*n,1));
Lower = Sabc(m*n*alpha);
Upper = Sabc(m*n*(1-alpha));

% find the clusters above the thresholds
XY = [X(find(abc<Lower | abc>Upper)) Y(find(abc<Lower | abc>Upper))];
XYZ = [XY ones(length(XY),1)];
A = spm_clusters(XYZ');
NCluster = max(A);
ClusterSize = zeros(NCluster,1);
for i = 1:NCluster
    ClusterSize(i) = length(find(A == i));
end
max(ClusterSize)
% find the number of clusters at different sizes and different thresholds
%%

figure(1)
hist(a,20);
figure(2)
hist(ab,20)
figure(3)
hist(abc,20)



        zh0 = norminv(length(find(bstat(:,j,k) < pointEst(j,k)))/nboot);
        zA = norminv(alpha/2);
        z1mA = norminv(1 - alpha/2);
        ThetaDiff = (sum(theta(:,j))/N) - theta(:,j);
        acc = (sum(ThetaDiff.^3))/(6*(sum(ThetaDiff.^2))^(3/2));
        Alpha1(j,k) = normcdf(zh0 + (zh0+zA)/(1 - acc*(zh0 + zA)));
        Alpha2(j,k) = normcdf(zh0 + (zh0+z1mA)/(1 - acc*(zh0 + z1mA)));

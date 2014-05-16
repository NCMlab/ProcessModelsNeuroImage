% THe idea is to create random fields, then smooth them.
m = 100;
n = 100;
p = 100;
[X Y Z]=meshgrid(1:n,1:m,1:p);
a = randn(m,n,p);

b = randn(m,n,p);

c = randn(m,n,p);
% smooth images
Sa = zeros(size(a));
Sb = zeros(size(b));
Sc = zeros(size(c));
spm_smooth(a,Sa,[8 8 8]);
spm_smooth(b,Sb,[8 8 8]);
spm_smooth(c,Sc,[8 8 8]);

alpha = 0.005;
% multiply them

abc = Sa.*Sb.*Sc;
% threshold them
Sabc = sort(reshape(abc,m*n*p,1));
Lower = Sabc(m*n*p*alpha);
Upper = Sabc(m*n*p*(1-alpha));

% find the clusters above the thresholds
XYZ = [X(find(abc<Lower | abc>Upper)) Y(find(abc<Lower | abc>Upper)) Z(find(abc<Lower | abc>Upper))];
A = spm_clusters(XYZ');
NCluster = max(A);
ClusterSize = zeros(NCluster,1);
for i = 1:NCluster
    ClusterSize(i) = length(find(A == i));
end
max(ClusterSize)
AllClusterSizes = unique(ClusterSize);
NumClusters = zeros(length(AllClusterSizes),1);
for i = 1:length(AllClusterSizes)
    NumClusters(i) = length(find(ClusterSize == AllClusterSizes(i)));
end
[AllClusterSizes NumClusters]    
    
% find the number of clusters at different sizes and different thresholds
%%

figure(1)
hist(a,20);
figure(2)
hist(ab,20)
figure(3)
hist(abc,20)



     
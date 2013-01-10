P = spm_select(Inf,'image');
for j = 1:size(P,1)
    V = spm_vol(deblank(P(j,:)));
    I = spm_read_vols(V);
    F = find(I);
    [x y z] = ind2sub(V.dim,F);
    L = [x y z];
    A = spm_clusters(L');
    Ncluster = max(A);
    ClusterSizes = zeros(Ncluster,1);
    ExtentThreshold = 50;
    thrString = sprintf('_Kthr%d',ExtentThreshold);
    for i = 1:Ncluster
        ClusterSizes(i) = length(find(A == i));
    end
    
    F2 = find(ClusterSizes > ExtentThreshold);
    NsignCluster = length(F2);
    threshI = zeros(size(I));
    for i = 1:NsignCluster
        threshI(F(find(A == F2(i)))) = 1;
    end
    Vo = V;
    [PathName FileName Ext] = fileparts(V.fname);
    Vo.fname = fullfile(PathName,[FileName thrString Ext]);
    spm_write_vol(Vo,threshI);
end
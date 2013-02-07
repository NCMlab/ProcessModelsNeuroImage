function OutImage = subfnApplyThresholdsToImages(InputImage,HeightThreshold,ExtentThreshold)
% P = spm_select(Inf,'image');
% HeightThreshold = 0.5;
% ExtentThreshold = 300;
% convert HeightThreshold value to a string

Hthr = num2str(HeightThreshold);
Hthr(findstr(Hthr,'.')) = 'p';


for j = 1:size(InputImage,1)
    % check to see if this image is a binary significance image
    if ~isempty(strfind(deblank(InputImage(j,:)),'sign0'))
        thrString = sprintf('_Kthr%d',ExtentThreshold);
    else
        thrString = sprintf('_Hthr%s_Kthr%d',Hthr,ExtentThreshold);
    end
    % load the image
    V = spm_vol(deblank(InputImage(j,:)));
    I = spm_read_vols(V);
    
    % find voxels above threhsold in the POSITIVE direction    
    F = find(I > HeightThreshold);
    % convert supreathreshold voxels to locations    
    [x y z] = ind2sub(V.dim,F);
    L = [x y z];
    % identify clusters    
    A = spm_clusters(L');
    Ncluster = max(A);
    ClusterSizes = zeros(Ncluster,1);
    % find the cluster sizes    
    for i = 1:Ncluster
        ClusterSizes(i) = length(find(A == i));
    end
    % Find the clusters that are "big enough"   
    F2 = find(ClusterSizes > ExtentThreshold);
    NsignCluster = length(F2);
    % create an image of suprathreshold clusters    
    threshI = zeros(size(I));
    for i = 1:NsignCluster
        threshI(F(find(A == F2(i)))) = 1;
    end
    % convert the suprathreshold voxel values to original values
    POSthresh = I.*threshI;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find voxels above threhsold in the NEGATIVE direction    
    F = find(I < -1*HeightThreshold);
    % convert supreathreshold voxels to locations    
    [x y z] = ind2sub(V.dim,F);
    L = [x y z];
    % identify clusters    
    A = spm_clusters(L');
    Ncluster = max(A);
    ClusterSizes = zeros(Ncluster,1);
    % find the cluster sizes    
    for i = 1:Ncluster
        ClusterSizes(i) = length(find(A == i));
    end
    % Find the clusters that are "big enough"   
    F2 = find(ClusterSizes > ExtentThreshold);
    NsignCluster = length(F2);
    % create an image of suprathreshold clusters    
    threshI = zeros(size(I));
    for i = 1:NsignCluster
        threshI(F(find(A == F2(i)))) = 1;
    end
    % convert the suprathreshold voxel values to original values
    NEGthresh = I.*threshI;
    OutI = POSthresh + NEGthresh;
    Vo = V;
    [PathName FileName Ext] = fileparts(V.fname);
    Vo.fname = fullfile(PathName,[FileName thrString Ext]);
    spm_write_vol(Vo,OutI);
    fprintf(1,'Image written to: %s\n,',Vo.fname);
end
OutImage = Vo.fname;
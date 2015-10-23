F = dir('Path_count*');
N = length(F);
% Load the first one to find sizes
load(F(1).name)
[m n o] = size(PermResults.beta);
NPaths = size(PermResults.Paths,1);

MaxTFCE_t = zeros(m,n,N);
MinTFCE_t = zeros(m,n,N);
MaxTFCE_PathsTNorm = zeros(NPaths,N);
MaxTFCE_PathsTNorm = zeros(NPaths,N);

load ../data/ModelInfo
Vo = ModelInfo.DataHeader;
Vo.fname = fullfile(pwd,'tempP.nii');
Vm = Vo;
Vm.fname = fullfile(pwd,'tempMask.nii');
I = zeros(Vo.dim);
I(ModelInfo.Indices) = 1;
spm_write_vol(Vm,I);
TFCEparams = ModelInfo.TFCEparams;
Nvoxels = o;

for k = 1:N
    % load up each permutation
    t = tic;
    load(F(k).name)
    for i = 1:m
        for j = 1:n
            if abs(PermResults.beta(i,j,1)) > 0

                I = zeros(Vo.dim);
                I(ModelInfo.Indices) = squeeze(PermResults.t(i,j,:));
                spm_write_vol(Vo,I);

                [MaxTFCE_t(i,j,k), MinTFCE_t(i,j,k)] = subfnApplyTFCE(Vo.fname,Vm.fname, TFCEparams);
            end
        end
    end
    % Now recalculate the paths
    for mm = 1:NPaths
        PathHolder = zeros(Nvoxels,1);
        for nn = 1:Nvoxels
            PathHolder(nn,1) = PermResults.PathsTnorm{nn}{mm};
        end
        I = zeros(Vo.dim);
        I(ModelInfo.Indices) = PathHolder(:,mm);
        spm_write_vol(Vo,I);
        [MaxTFCE_PathsTNorm(mm,k), MinTFCE_PathsTNorm(mm,k)] = subfnApplyTFCE(Vo.fname,Vm.fname, TFCEparams);
    end
    fprintf(1,'Recalculating perm %d of %d in %0.2f sec\n',k,N,toc(t));
end

save RecalutePermTFCEv2 MaxTFCE_PathsTNorm MinTFCE_PathsTNorm MaxTFCE_t MinTFCE_t
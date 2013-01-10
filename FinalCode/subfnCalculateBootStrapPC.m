function subfnCalculateBootStrapPC(InDataFile)
load(InDataFile)
% reset the random number seed
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% Find out what count this is based on the job ID
[a b] = unix('echo $SGE_TASK_ID');
count = str2num(b);
fprintf(1,'Working on Block %d\n',count);


NSub = length(data.X);
remove_row_means = 1;
BootStrapResamples = zeros(NSub,Nboot,'uint16');
if ~isempty(data.STRAT)
    Gr1 = find(data.STRAT == 0);
    NGr1 = length(Gr1);
    Gr2 = find(data.STRAT == 1);
    NGr2 = length(Gr2);
else
    NGr1 = [];
    NGr2 = [];
end
for i = 1:Nboot
    if isempty(data.STRAT)
        Samp =  floor(N*rand(N,1))+1;
    else
        Samp1 = floor(NGr1*rand(NGr1,1))+1;
        Samp2 = floor(NGr2*rand(NGr2,1))+1+NGr1;
        Samp = [Samp1; Samp2];
    end
    BootStrapResamples(:,i) = Samp;
end
NVox = size(data.M,3);
BootStrapResampleImages = zeros(NVox,Nboot);
for i = 1:Nboot
    TrainData.X = data.X(BootStrapResamples(:,i));
    TrainData.M = data.M(BootStrapResamples(:,i),:,:);
    TrainData.Y = data.Y(BootStrapResamples(:,i));
    TrainData.ModelNum = data.ModelNum;
    % Apply PCA
    [lambdas, eigenimages_noZeroes, w] = pca_f(squeeze(TrainData.M)', remove_row_means);
    % calculate the subject scaling factors
    ssf = squeeze(TrainData.M) * eigenimages_noZeroes;
    ssfSubset = ssf;
    behav_fit_coef = subfnregress(TrainData.Y, [ssfSubset(:,selected_PCs) TrainData.X]);
    % create the SSF image
    temp = eigenimages_noZeroes(:, selected_PCs) * behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent
    % Created the weighted sum pf SSFs regressor
    M = ssfSubset(:,selected_PCs)*behav_fit_coef(1 + 1:1 + length(selected_PCs));
    a_coef = subfnregress(M,TrainData.X);
    % Calculate the mediated effect
    temp = temp.*a_coef(end);
    behav_fit_composite_PC_image = temp / norm(temp);
    clear temp;
    BootStrapResampleImages(:,i) = behav_fit_composite_PC_image;
      if ~mod(i,10)
        fprintf(1,'Bootstrap resample: %5d\n',i);
    end
end

outFile = fullfile(OutDir,['BootStrapChunk' num2str(Nboot) '_' num2str(count)])
Str = sprintf('save %s %s' , outFile, 'BootStrapResampleImages');
eval(Str);
fprintf(1,'Saved file: %s\n',outFile);
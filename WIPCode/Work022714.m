function WIPLoadBootStraps
%%
if nargin == 0
    SelectedDir = spm_select(1,'dir');
end
cd(SelectedDir)
cd('BootStraps');



F = dir('BootStrapResults_*.mat');
NJobs = length(F);



clear BootStrap
i = 1;
tic
load(F(i).name);
NPaths = size(BootStrap.Paths(:,1,1),1);
NVox = size(BootStrap.Paths(1,:,1),2);
NBootSplit = size(BootStrap.Paths(1,1,:),3);
Nboot = NBootSplit*NJobs;

MaxPaths = zeros(Nboot,NPaths);
MinPaths = zeros(Nboot,NPaths);



i = 1;
Paths = cell2mat(BootStrap.Paths);
for j = 1:NPaths
    for k = 1:NBootSplit
        MaxPaths((i-1)*NBootSplit+1:(i-1)*NBootSplit+NBootSplit,j) = max(squeeze(Paths(j,:,k)));
        MinPaths((i-1)*NBootSplit+1:(i-1)*NBootSplit+NBootSplit,j) = min(squeeze(Paths(j,:,k)));
    end
end

t = toc;
fprintf(1,'Finished loading job %d of %d in %0.1f seconds.\n',i,NJobs,t);
for i = 2:NJobs
    tic
    clear BootStrap
    load(F(i).name)
    Paths = cell2mat(BootStrap.Paths);
    for j = 1:NPaths
        for k = 1:NBootSplit
            MaxPaths((i-1)*NBootSplit+1:(i-1)*NBootSplit+NBootSplit,j) = max(squeeze(Paths(j,:,k)));
            MinPaths((i-1)*NBootSplit+1:(i-1)*NBootSplit+NBootSplit,j) = min(squeeze(Paths(j,:,k)));
        end
    end
    t = toc;
    fprintf(1,'Finished loading job %d of %d in %0.1f seconds.\n',i,NJobs,t);
end
clear BootStrap

%%
% load point estimates
PEPaths = zeros(NPaths,NVox);
for i = 1:NVox
    for j = 1:NPaths
        PEPaths(j,i) = PointEstResults{i}.Paths{j};
    end
end


BSPaths = zeros(NPaths,NVox);
for i = 1:NVox
    for j = 1:NPaths
        BSPaths(j,i) = BootParameters{i}.Paths{j};
    end
end


a = 0.05;
c = floor(a*Nboot);
Limits = zeros(NPaths,2);
for j = 1:NPaths
    S = sort(MaxPaths(:,j),'descend');
    Limits(j,1) = S(c);
    S = sort(MinPaths(:,j),'ascend');
    Limits(j,2) = S(c);
end
for j = 1:NPaths;
    figure(j)
    clf
    hist(PEPaths(j,:),50);
    h = line([Limits(j,1) Limits(j,1)],[0 1000]);
    set(h,'Color','r')
    h = line([Limits(j,2) Limits(j,2)],[0 1000]);
    set(h,'Color','r')
end

%%
Nvox = 10001;
PEa = zeros(Nvox,1);
PEb = zeros(Nvox,1);
PEc = zeros(Nvox,1);
for i = 1:Nvox
    PEa(i) = PointEstResults{i}.beta(2,2);
    PEb(i) = PointEstResults{i}.beta(3,3);
    PEc(i) = PointEstResults{i}.beta(4,4);
end
PEPath = PEa.*PEb.*PEc;


a = squeeze(BETAS(2,2,:,:));
b = squeeze(BETAS(3,3,:,:));
c = squeeze(BETAS(4,4,:,:));

length(unique(a(100,:)))
length(unique(b(100,:)))
length(unique(c(100,:)))

PATHS = a.*b.*c;
maxPATHS = max(PATHS);
SmaxPATHS = sort(maxPATHS,'descend');

minPATHS = min(PATHS);
SminPATHS = sort(minPATHS,'ascend');


alp = 0.025;
c = alp*Nboot;
Limit1 = SmaxPATHS(c)
Limit2 = SminPATHS(c)


figure(1)
clf
hist(maxPATHS,50)
h = line([Limit1 Limit1],[0 50]);
set(h,'Color','r')
figure(2)
clf
hist(minPATHS,50)
h = line([Limit2 Limit2],[0 50]);
set(h,'Color','r')

length(find(PEPath>Limit1 | PEPath<Limit2));

figure(3)
hist(PEPath,50)
h = line([Limit2 Limit2],[0 50]);
set(h,'Color','r')
h = line([Limit1 Limit1],[0 50]);
set(h,'Color','r')

%%
% Find Max
P1 = squeeze(Paths(1,:,:));
P2 = squeeze(Paths(1,:,:));
P3 = squeeze(Paths(1,:,:));
MaxPath1 = max(P1);
MaxPath2 = max(P2);
MaxPath3 = max(P3);
figure
%hist(MaxPath1,20)
hist(MaxPath1,20)
hist(MaxPath3,20)
%%
a = squeeze(BootStrap.beta(2,2,:,1));
b = squeeze(BootStrap.beta(3,3,:,1));
c = squeeze(BootStrap.beta(4,4,:,1));

P1 = a.*b.*c;

% 
a = squeeze(BETAS(2,2,:,1:15*5));
b = squeeze(BETAS(3,3,:,1:15*5));
c = squeeze(BETAS(4,4,:,1:15*5));
P1 = a.*b.*c;

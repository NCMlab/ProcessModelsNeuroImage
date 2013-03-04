function [SPM xSPM] = DisplayCov(varargin)
%clear classes
if nargin == 3
    HeightThr = varargin{2};
    ExtentThr = varargin{3};
    P = varargin{1};
elseif nargin == 1
    P = varargin{1};
    HeightThr = [];
    ExtentThr = [];
elseif nargin < 1
    % Select image
    P = spm_select(1,'image','Select Image')
    HeightThr = [];
    ExtentThr = [];
end

% Filename DisplayCov.m
% Written by: Jason Steffener
% Date: 09/07
%
%-----------------------------------------------------------------------
u      = -Inf;
k      = 0;
n = 1;
STAT = 'T';
df = [1 100]; % Need degreees of freedom
V = spm_vol(P);
Data = spm_read_vols(spm_vol(P));
%
% Select direction of inference and calculate the max and min values
MinRange = round(min(min(min(Data)))*100)/100;
MaxRange = round(max(max(max(Data)))*100)/100;
if isempty(HeightThr) 
    DirectionStr = ['Display:' '[range:' num2str(MinRange) ':' num2str(MaxRange) ']']
    str = 'positive|negative';
    direction = spm_input(DirectionStr,'+1','b',str,[],1);
    u  = spm_input(['threshold {',STAT,' or p value}'],'+1','r',3.09,1);
else
    if sign(HeightThr) > 0
        direction = 'positive';
    else 
        direction = 'negative';
    end
    u = HeightThr;
end

%
% do some error checking to make sure values are within range
if ~isempty(strmatch(direction,'negative'))
    while ~((u<=0)&(u>MinRange))
        warndlg('Please enter a value within range','','modal')
        u  = spm_input(['threshold {',STAT,' or p value}'],'+1','r',3.09,1);
    end
else
    while ~((u>=0)&(u<MaxRange))
        warndlg('Please enter a value within range','','modal')
        u  = spm_input(['threshold {',STAT,' or p value}'],'+1','r',3.09,1);
    end
end

thresDesc = [STAT '=' num2str(u) ];

DIM = V.dim;
[X Y Z] = ndgrid([1:V.dim(1)],[1:V.dim(2)],[1:V.dim(3)]);

% Check direction and if negative flip image around before thresholding
if strmatch(direction,'negative')
    Data = Data.*-1;
end
F  = find(Data > abs(u));
ZValues = zeros(1,length(F));
%

%
XYZ = zeros(3,length(F));
XYZmm = zeros(3,length(F));
SliceSize = DIM(1)*DIM(2);
for i = 1:length(F)
    XYZ(:,i) = [rem(rem(F(i),SliceSize),DIM(1)) ceil(rem(F(i),SliceSize)/DIM(1)) ceil(F(i)/SliceSize)];
     ZValues(1,i) = Data(XYZ(1,i), XYZ(2,i), XYZ(3,i));
end
% Select cluster size threshold based on number of pixels
if isempty(ExtentThr)
    k     = spm_input('& extent threshold {voxels}','+1','r',0,1,[0,Inf]);
else
    k = ExtentThr;
end
    
A     = spm_clusters(XYZ);
Q     = [];
% This selects cluster sizes. It places a threshold on the cluster size to
% keep. Any clusters of a size lower than k are excluded.
for i = 1:max(A)
    j = find(A == i);
    if length(j) >= k; 
        Q = [Q j]; 
    end
end
ZValues     = ZValues(:,Q);


XYZ   = XYZ(:,Q);
%      temp = [XYZ(:,i); 1]'*V.mat';
%      XYZmm(:,i) = temp(1:3)';

XYZmm = V.mat(1:3,:)*[XYZ; ones(1,size(XYZ,2))];
if isempty(Q)
    warning(sprintf('No voxels survive extent threshold k=%0.2g',k))
end
    units = {'mm' 'mm' 'mm'};
% [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');
% hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
% hMIPax = spm_mip_ui(ZValues,XYZmm,V.mat,V.dim',hMIPax,units);
[FWHM,VRpv] = spm_est_smoothness(P,P);
% This is a special catch 
if sum(FWHM == [Inf Inf Inf])==3
    FWHM = [20 20 20];
end
%FWHM = [8 8 8];
R           = spm_resels_vol(spm_vol(P),FWHM)';
%R = [1 1 1 1];
M = V.mat;
DIM = V.dim';
%    hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,xSPM.M,xSPM.DIM,hMIPax,xSPM.units);
xSPM = {};
xSPM.swd = pwd;%fileparts(P);
xSPM.Z = ZValues;
xSPM.n = 1;
xSPM.STAT = STAT;
xSPM.df = df;
xSPM.u = u;
xSPM.k = k;
xSPM.XYZ = XYZ;
xSPM.XYZmm = XYZmm;
xSPM.S = length(XYZ);
xSPM.R = R;
xSPM.FWHM = FWHM;
xSPM.M = M;
xSPM.VOX = diag(V.mat(1:3,1:3))';
xSPM.Vspm = V;
xSPM.Ps = [];
xSPM.thresDesc = thresDesc;
xSPM.STATstr = thresDesc;
% check direction and change the title accordingly
if strmatch(direction,'negative')
    xSPM.title = 'Negative Direction';
else
    xSPM.title = 'Positive Direction';
end
xSPM.Ic = 1;
xSPM.uc = [];
xSPM.Pc = [];
xSPM.Pp = [];
xSPM.Im = [];
xSPM.pm = [];
xSPM.Ex = [];
xSPM.VRpv = 0;
xSPM.DIM = DIM;
%%% need to build up the contrast matrix
SPM = {};
SPM.xVol = {};
SPM.xVol.M = M;
SPM.xVol.iM = inv(M);
SPM.xVol.DIM = DIM;
SPM.xVol.units = units;
SPM.swd = pwd;%fileparts(P);
SPM.xX.X = [1];
SPM.xX.nKX = SPM.xX.X;
SPM.xX.name = 'Cov';
SPM.xCon = [];
SPM.xY.VY = V;
SPM.xY.P = P;
SPM.xX.xKXs.X = SPM.xX.X;
SPM.xX.xKXs.tol= 2.8329e-13;
SPM.xX.xKXs.ds= 1;
SPM.xX.xKXs.u= 1;
SPM.xX.xKXs.v= 1;
SPM.xX.xKXs.rk= 1;
SPM.xX.xKXs.oP= [];
SPM.xX.xKXs.oPp= [];
SPM.xX.xKXs.ups= [];
SPM.xX.xKXs.sus= [];
SPM.xCon = {};
SPM.xCon(1).name = 'Contrast 1';
SPM.xCon(1).STAT = STAT;
SPM.xCon(1).c = 1;
SPM.xCon(1).X0 = {};
SPM.xCon(1).iX0 = 'c';
SPM.xCon(1).X1o = [];
SPM.xCon(1).eidf =1;
SPM.xCon(1).Vcon = V;
SPM.xCon(1).Vspm = V;
% Required fields of SPM
%
% xVol   - structure containing details of volume analysed
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xCon   - Contrast definitions structure array
%        - (see also spm_FcUtil.m for structure, rules & handling)
% .name  - Contrast name
% .STAT  - Statistic indicator character ('T', 'F' or 'P')
% .c     - Contrast weights (column vector contrasts)
% .X0    - Reduced design matrix data (spans design space under Ho)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          (Matrix in SPM99b)  Extract using X0 = spm_FcUtil('X0',...
% .iX0   - Indicates how contrast was specified:
%          If by columns for reduced design matrix then iX0 contains the
%          column indices. Otherwise, it's a string containing the
%          spm_FcUtil 'Set' action: Usually one of {'c','c+','X0'}
% .X1o   - Remaining design space data (X1o is orthogonal to X0)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          (Matrix in SPM99b)  Extract using X1o = spm_FcUtil('X1o',...
% .eidf  - Effective interest degrees of freedom (numerator df)
%        - Or effect-size threshold for Posterior probability
% .Vcon  - Name of contrast (for 'T's) or ESS (for 'F's) image
% .Vspm  - Name of SPM image

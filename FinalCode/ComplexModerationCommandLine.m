function [Model1, Results] = ComplexModerationCommandLine(A,B,C,D,COV,names,alpha,Nboot,Model)
switch nargin
    case 0
        error('MediationCommandLine needs at least 3 input arguments')
    case 1
        error('MediationCommandLine needs at least 3 input arguments')
    case 2
        error('MediationCommandLine needs at least 3 input arguments')
    case 3
        D = [];
        COV = [];
        names = {'A' 'B' 'C'};
        alpha = 0.05;
        Nboot = 1000;
        ModelFlag = 0;
    case 4
        COV = [];
        names = {'A' 'B' 'C' 'D'};
        alpha = 0.05;
        Nboot = 1000;
        ModelFlag = 0;
    case 5
        names = {'A' 'B' 'C' 'D'};
        if length(names) < size(COV,2) + 4
            for i = 1:size(COV,2)
                names{4+i} = sprintf('COV%03d',i);
            end
        end
        alpha = 0.05;
        Nboot = 1000;
        ModelFlag = 0;
    case 6
        if isempty(names)
            names = {'A' 'B' 'C' 'D'};
        end
        if length(names) < size(COV,2) + 4
            for i = 1:size(COV,2)
                names{4+i} = sprintf('COV%03d',i);
            end
        end
        alpha = 0.05;
        Nboot = 1000;
        ModelFlag = 0;
    case 7
        if isempty(names)
            names = {'A' 'B' 'C' 'D'};
        end
        if length(names) < size(COV,2) + 4
            for i = 1:size(COV,2)
                names{4+i} = sprintf('COV%03d',i);
            end
        end
        Nboot = 1000;
        ModelFlag = 0;
    case 8
        if isempty(names)
            names = {'A' 'B' 'C' 'D'};
        end
        if length(names) < size(COV,2) + 4
            for i = 1:size(COV,2)
                names{4+i} = sprintf('COV%03d',i);
            end
        end
        ModelFlag = 0;
    case 9
        if isempty(names)
            names = {'A' 'B' 'C' 'D'};
        end
        if length(names) < size(COV,2) + 4
            for i = 1:size(COV,2)
                names{4+i} = sprintf('COV%03d',i);
            end
        end
        ModelFlag = 1;
end

data = [A B C D COV];
[NSub Nvar] = size(data);
NCov = size(COV,2);

Model1 = {};
Model1.BaseDir = pwd;
Model1.Names = names;
Model1.data = data;
Model1.Nboot = Nboot;
Model1.Nperm = 0;

Model1.Indices = 1;
Model1.NJobSplit = 1;
Model1.Thresholds = alpha;
% Startification is used when the resamples are created and needs to be a
% binomial parameter for right now. The use of a stratification variable is
% so that when the reampling is performed each resample maintains the
% number of subjects as in the stratification parameter. This is most
% applicable when there are multiple groups with different sample sizes.
Model1.STRAT = [];
Model1.Nsub = NSub;
Model1.Nvar = Nvar;
Model1.Nvoxels = 1;

% Prepare the output data header
DataHeader.fname = '';
DataHeader.descrip = '';
DataHeader.dt = [16 0];
Model1.DataHeader = DataHeader;
if ~ModelFlag
    % Create the direct effects model assuming the covariates are for all
    % variables
    Direct = zeros(Nvar);
    Direct(1,[3]) = 1;
    Direct(2,[3]) = 1;
    for i = 1:NCov
        Direct(3+i,[3]) = 1;
    end
    Inter = zeros(Nvar);
    Inter([1 2],3) = 1;
    Paths = zeros(Nvar);
    %Paths(1,2) = 1;
    Paths(1,3) = 1;
else
    Direct = Model.Direct;
    Inter = Model.Inter;
    Paths = Model.Paths;
end
Model1.Direct = Direct;
Model1.Inter = Inter;
Model1.Paths = Paths;
Results = OneVoxelProcessBootstrap(Model1);
PrintResults(Model1,Results)


function subfnPrintFSResults_beta(SelectedPath,alpha,fid,StandardizeFlag,FDRflag)
% write Fressurfer results to tables
switch nargin
    case 3
        StandardizeFlag = 0;
        FDRflag = 0;
    case 2
        fid = 1;
        StandardizeFlag = 0;
        FDRflag = 0;
    case 1
        alpha = 0.05;
        fid = 1;
        StandardizeFlag = 0;
        FDRflag = 0;
    case 0
        SelectedPath = spm_select(1,'dir');
        alpha = 0.05;
        fid = 1;
        StandardizeFlag = 0;
        FDRflag = 0;
end

cd(SelectedPath)
if exist('AnalysisParameters.mat')
    load AnalysisParameters
    
    F = dir('Results_*.mat');
    if length(F) == AnalysisParameters.NJobSplit
        NvoxelsPerJob = ceil(AnalysisParameters.Nvoxels/AnalysisParameters.NJobSplit);
        % prealocate memory for AllParameters
        AllParameters = cell(AnalysisParameters.Nvoxels,1);
        % load up all the data
        for i = 1:length(F)-1
            clear Parameters
            load(F(i).name)
            AllParameters((i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob) = Parameters;
        end
        clear Parameters
        load(F(i+1).name)
        AllParameters(i*NvoxelsPerJob+1:end) = Parameters;
        if StandardizeFlag
            
            Fdata = dir('data_*.mat');
            AllData = cell(AnalysisParameters.Nvoxels,1);
            % load up all the data
            for i = 1:length(F)-1
                clear data
                load(Fdata(i).name)
                AllData{(i-1)*NvoxelsPerJob + 1:i*NvoxelsPerJob} = data;
            end
            clear data
            load(Fdata(i+1).name)
            AllData{i*NvoxelsPerJob+1:end} = data;
        end
    else
        errordlg('This analysis is not finished.');
    end
else
    errordlg('this folder does not have the required AnalysisParameters.mat file');
end
[PathName FileName] = fileparts(SelectedPath);
% Setup FS names
[outNames Hemi MeasurementType] = subfnModifyFreesurferNames(AnalysisParameters.Header);
%%
Nvoxels = AnalysisParameters.Nvoxels;
% for i = 1:Nvoxels
%     fprintf(1,'%d\t%s\n',i,AnalysisParameters.Header{i});
% end
%

ModelNum = AnalysisParameters.ModelNum;
thr = 1.96;
%alpha = 0.05;%/Nvoxels;
% fin dthe closest threshold to the calculated alpha
D = AnalysisParameters.Thresholds - alpha;
TempAlpha = AnalysisParameters.Thresholds(find(D==min(abs(D))));
TempAlpha = TempAlpha(1);
StrAlpha = num2str(TempAlpha);
Dot=findstr(StrAlpha,'.');
AlphaLimit = ['alpha' StrAlpha(Dot+1:end)];

fprintf(fid,'%s\n',FileName);


PrintHeaderFlag = 1;
% %% Find FDR rates
% switch ModelNum
%     case '58'
%         % Get p values
%         a = zeros(Nvoxels,1);
%         p  = zeros(Nvoxels,1);
%         w = zeros(Nvoxels,1);
%         cP = zeros(Nvoxels,1);
%         b = zeros(Nvoxels,1);
%         q = zeros(Nvoxels,1);
%         v = zeros(Nvoxels,1);
%         for i = 1:Nvoxels
%             P = AllParameters{i};
%             % Model 1
%             Sa = getfield(P.Model1{1},[P.names.X]);
%             a(i) = Sa.p;
%             Sp = getfield(P.Model1{1},[P.names.W]);
%             p(i) = Sp.p;
%             Sw = getfield(P.Model1{1},[P.names.X '_x_' P.names.W]);
%             w(i) = Sw.p;
%             % Model 2
%             ScP =getfield(P.Model2,[P.names.X]);
%             cP(i) = ScP.p;
%             Sb = getfield(P.Model2,[P.names.M{1}]);
%             b(i) = Sb.p;
%             Sq = getfield(P.Model2,[P.names.W]);
%             q(i) = Sq.p;
%             Sv = getfield(P.Model2,[P.names.M{1} '_x_' P.names.W]);
%             v(i) = Sv.p;
%         end
%         FDRa = FDR(a,alpha,FDRflag);
%         FDRp = FDR(p,alpha,FDRflag);
%         FDRw = FDR(w,alpha,FDRflag);
%         FDRcP = FDR(cP,alpha,FDRflag);
%         FDRb = FDR(b,alpha,FDRflag);
%         FDRq = FDR(q,alpha,FDRflag);
%         FDRv = FDR(v,alpha,FDRflag);
%         if isempty(FDRa)
%             FDRa = 10^-10;
%         end
%         if isempty(FDRp)
%             FDRp = 10^-10;
%         end
%         if isempty(FDRw)
%             FDRw = 10^-10;
%         end
%         if isempty(FDRcP)
%             FDRcP = 10^-10;
%         end
%         if isempty(FDRb)
%             FDRb = 10^-10;
%         end
%         if isempty(FDRq)
%             FDRq = 10^-10;
%         end
%         if isempty(FDRv)
%             FDRv = 10^-10;
%         end
%         if FDRflag
%             wTHR = FDRv;
%         else
%             %wTHR = alpha;
%             wTHR = 0.05;
%         end
%         if FDRflag
%             vTHR = FDRv;
%         else
%             vTHR = 0.05;
%             %vTHR = alpha;
%         end
%     case '14'
%         % Get p values
%         a = zeros(Nvoxels,1);
%         cP = zeros(Nvoxels,1);
%         b = zeros(Nvoxels,1);
%         q = zeros(Nvoxels,1);
%         v = zeros(Nvoxels,1);
%         for i = 1:Nvoxels
%             P = AllParameters{i};
%             % Model 1
%             Sa = getfield(P.Model1{1},[P.names.X]);
%             a(i) = Sa.p;
%             % Model 2
%             ScP =getfield(P.Model2,[P.names.X]);
%             cP(i) = ScP.p;
%             Sb = getfield(P.Model2,[P.names.M{1}]);
%             b(i) = Sb.p;
%             Sq = getfield(P.Model2,[P.names.V]);
%             q(i) = Sq.p;
%             Sv = getfield(P.Model2,[P.names.M{1} '_x_' P.names.V]);
%             v(i) = Sv.p;
%         end
%         FDRa = FDR(a,alpha,FDRflag);
%         FDRcP = FDR(cP,alpha,FDRflag);
%         FDRb = FDR(b,alpha,FDRflag);
%         FDRq = FDR(q,alpha,FDRflag);
%         FDRv = FDR(v,alpha,FDRflag);
%         if isempty(FDRa)
%             FDRa = 10^-10;
%         end
%         if isempty(FDRcP)
%             FDRcP = 10^-10;
%         end
%         if isempty(FDRb)
%             FDRb = 10^-10;
%         end
%         if isempty(FDRq)
%             FDRq = 10^-10;
%         end
%         if isempty(FDRv)
%             FDRv = 10^-10;
%         end
%         if FDRflag
%             vTHR = FDRv;
%         else
%             vTHR = 0.05;
%             %vTHR = alpha;
%         end
%     case '7'
%         % Get p values
%         a = zeros(Nvoxels,1);
%         p  = zeros(Nvoxels,1);
%         w = zeros(Nvoxels,1);
%         cP = zeros(Nvoxels,1);
%         b = zeros(Nvoxels,1);
%         for i = 1:Nvoxels
%             P = AllParameters{i};
%             % Model 1
%             Sa = getfield(P.Model1{1},[P.names.X]);
%             a(i) = Sa.p;
%             Sp = getfield(P.Model1{1},[P.names.W]);
%             p(i) = Sp.p;
%             Sw = getfield(P.Model1{1},[P.names.X '_x_' P.names.W]);
%             w(i) = Sw.p;
%             % Model 2
%             ScP =getfield(P.Model2,[P.names.X]);
%             cP(i) = ScP.p;
%             Sb = getfield(P.Model2,[P.names.M{1}]);
%             b(i) = Sb.p;
%         end
%         FDRa = FDR(a,alpha,FDRflag);
%         FDRp = FDR(p,alpha,FDRflag);
%         FDRw = FDR(w,alpha,FDRflag);
%         FDRcP = FDR(cP,alpha,FDRflag);
%         FDRb = FDR(b,alpha,FDRflag);
%         if isempty(FDRa)
%             FDRa = 10^-10;
%         end
%         if isempty(FDRp)
%             FDRp = 10^-10;
%         end
%         if isempty(FDRw)
%             FDRw = 10^-10;
%         end
%         if isempty(FDRcP)
%             FDRcP = 10^-10;
%         end
%         if isempty(FDRb)
%             FDRb = 10^-10;
%         end
%         if FDRflag
%             wTHR = FDRw;
%         else
%             wTHR = 0.05;
% %             wTHR = alpha;
%         end
%     case '4'
%         % Get p values
%         a = zeros(Nvoxels,1);
%         cP = zeros(Nvoxels,1);
%         b = zeros(Nvoxels,1);
%         for i = 1:Nvoxels
%             P = AllParameters{i};
%             % Model 1
%             Sa = getfield(P.Model1{1},[P.names.X]);
%             a(i) = Sa.p;
%             % Model 2
%             ScP =getfield(P.Model2,[P.names.X]);
%             cP(i) = ScP.p;
%             Sb = getfield(P.Model2,[P.names.M{1} '1']);
%             b(i) = Sb.p;
%         end
%         FDRa = FDR(a,alpha,FDRflag);
%         FDRcP = FDR(cP,alpha,FDRflag);
%         FDRb = FDR(b,alpha,FDRflag);
%         if isempty(FDRa)
%             FDRa = 10^-10;
%         end
%         if isempty(FDRcP)
%             FDRcP = 10^-10;
%         end
%         if isempty(FDRb)
%             FDRb = 10^-10;
%         end
% end
% 
% 

wTHR = 1;%0.05;
vTHR = 1;%0.05;
for i = 1:AnalysisParameters.Nvoxels
    P = AllParameters{i};
    %D = AllData{i};
    switch ModelNum
        case '58'
            NCondSteps = length(P.CondAB1);
            if PrintHeaderFlag
                fprintf(fid,'Interaction threshold = %0.4f, %0.4f\n',wTHR,vTHR);
                fprintf(fid,'%4s,%-30s,%-10s,%-10s,%7s,%7s,%7s,%7s,%7s,%7s,%7s,','ind','Region','Hemi','Measure','a','p','w','cP','b','q','v');
                ProbeValues = [0 10 25 34 50 66 75 90];
                for j = 2:NCondSteps
                    fprintf(fid,'%7.0f,',ProbeValues(j));
                end
                fprintf(fid,'\n');
                PrintHeaderFlag = 0;
            end
            % Check to make sure both interactions are significant
            % Model 1
            a = getfield(P.Model1{1},[P.names.X]);
            p = getfield(P.Model1{1},[P.names.W]);
            w = getfield(P.Model1{1},[P.names.X '_x_' P.names.W]);
            % Model 2
            cP =getfield(P.Model2,[P.names.X]);
            b = getfield(P.Model2,[P.names.M{1}]);
            q = getfield(P.Model2,[P.names.W]);
            v = getfield(P.Model2,[P.names.M{1} '_x_' P.names.W]);
            
            % Are BOTH interactions significant?
            if abs(w.p) < wTHR & abs(v.p) < vTHR
                % get the conditional effects
                CondEffectsSIGN = zeros(NCondSteps-1,1);
                CondEffects = zeros(NCondSteps-1,1);
                for j = 2:NCondSteps
                    Limits = getfield(P.CondAB1{j}.BCaci,AlphaLimit);
                    CondEffects(j-1) = P.CondAB1{j}.pointEst;
                    if prod(Limits) > 0
                        CondEffectsSIGN(j-1) = sign(Limits(1));
                        
                    end
                end
                if ~isempty(find(CondEffectsSIGN~=0))
                    % print out results
                    %fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,num2str(i),'?','area');
                    fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,outNames{i},Hemi{i},MeasurementType{i});
                    str = subfnCreateParamStr(a,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(p,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(w,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(cP,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(b,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(q,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(v,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    
                    for j = 2:NCondSteps
                        if CondEffectsSIGN(j-1) ~= 0
                            fprintf(fid,'%7.3f*,',CondEffects(j-1));
                        else
                            fprintf(fid,'%7.3f,',CondEffects(j-1));
                        end
                    end
                    fprintf(fid,'\n');
                end
            end
        case '14'
            NCondSteps = length(P.CondAB1);
            if PrintHeaderFlag
                fprintf(fid,'Interaction threshold = %0.4f\n',vTHR);
                fprintf(fid,'%4s,%-30s,%-10s,%-10s,%7s,%7s,%7s,%7s,%7s,','ind','Region','Hemi','Measure','a','cP','b','q','v');
                ProbeValues = [0 10 25 34 50 66 75 90];
                for j = 2:NCondSteps
                    fprintf(fid,'%7.0f,',ProbeValues(j));
                end
                fprintf(fid,'\n');
                PrintHeaderFlag = 0;
            end
            % Check to make sure both interactions are significant
            % Model 1
%             X = getfield(D,'X');
%             M = getfield(D,'M');
%             Y = getfield(D,'Y');
%             V = getfield(D,'V');
%             std_X = std(X);
%             std_M = std(M);
%             std_Y = std(Y);
%             std_V = std(V);
%             std_Inter = std(M.*V);
            a = getfield(P.Model1{1},[P.names.X]);
%             Sa = a.beta*std_X/std_M;
            % Model 2
            cP =getfield(P.Model2,[P.names.X]);
%             ScP = cP.beta*std_X/std_Y;
            b = getfield(P.Model2,[P.names.M{1}]);
%             Sb = b.beta*std_M/std_Y;
            q = getfield(P.Model2,[P.names.V]);
%             Sq = q.beta*std_V/std_Y;
            v = getfield(P.Model2,[P.names.M{1} '_x_' P.names.V]);
%             Sv = v.beta*std_Inter/std_Y;
            % Is the interaction significant?
            if abs(v.p) < vTHR
                % get the conditional effects
                CondEffectsSIGN = zeros(NCondSteps-1,1);
                %CondEffects = zeros(NCondSteps-1,1);
                CondEffects = {};
                %CondEffects = zeros(NCondSteps-1,1);
                for j = 2:NCondSteps
                    Limits = getfield(P.CondAB1{j}.BCaci,AlphaLimit);
                    %CondEffects(j-1) = P.CondAB1{j}.pointEst;
                    CondEffects{j-1}.beta = P.CondAB1{j}.pointEst;
                    CondEffects{j-1}.se = P.CondAB1{j}.bootSE;
                    if prod(Limits) > 0
                        CondEffectsSIGN(j-1) = sign(Limits(1));
                        CondEffects{j-1}.p = 0;
                    else
                        CondEffects{j-1}.p = 1;
                    end
                end
                if ~isempty(find(CondEffectsSIGN~=0))
                    % print out results
                    %fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,num2str(i),'?','area');
                    fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,outNames{i},Hemi{i},MeasurementType{i});
                    str = subfnCreateParamStr(a,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(cP,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(b,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(q,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(v,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    for j = 2:NCondSteps
                        %if CondEffectsSIGN(j-1) ~= 0
                            str = subfnCreateParamStr(CondEffects{j-1},alpha,FDRflag);
                            fprintf(fid,'%s,',str);
                        %else
%                             fprintf(fid,'%7.3f*,',CondEffects(j-1));
%                         else
%                             fprintf(fid,'%7.3f,',CondEffects(j-1));
%                         end
                        %end
                    end
                    fprintf(fid,'\n');
                end
            end
        case '7'
            NCondSteps = length(P.CondAB1);
            if PrintHeaderFlag
                fprintf(fid,'Interaction threshold = %0.4f\n',wTHR);
                fprintf(fid,'%4s,%-30s,%-10s,%-10s,%7s,%7s,%7s,%7s,%7s,','ind','Region','Hemi','Measure','a','p','w','cP','b');
                ProbeValues = [0 10 25 34 50 66 75 90];
                for j = 2:NCondSteps
                    fprintf(fid,'%7.0f,',ProbeValues(j));
                end
                fprintf(fid,'\n');             
                PrintHeaderFlag = 0;
            end
            % Check to make sure both interactions are significant
            % Model 1
            a = getfield(P.Model1{1},[P.names.X]);
            p = getfield(P.Model1{1},[P.names.W]);
            w = getfield(P.Model1{1},[P.names.X '_x_' P.names.W]);
            % Model 2
            cP =getfield(P.Model2,[P.names.X]);
            b = getfield(P.Model2,[P.names.M{1}]);
            % Is the interaction significant?
            if abs(w.p) < wTHR
                % get the conditional effects
                CondEffectsSIGN = zeros(NCondSteps-1,1);
                CondEffects = zeros(NCondSteps-1,1);
                for j = 2:NCondSteps
                    Limits = getfield(P.CondAB1{j}.BCaci,AlphaLimit);
                    CondEffects(j-1) = P.CondAB1{j}.pointEst;
                    if prod(Limits) > 0
                        CondEffectsSIGN(j-1) = sign(Limits(1));
                        
                    end
                end
                
                if ~isempty(find(CondEffectsSIGN~=0))
                    % print out results
                    %fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,num2str(i),'?','area');
                    fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,outNames{i},Hemi{i},MeasurementType{i});
                    str = subfnCreateParamStr(a,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(p,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(w,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(cP,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(b,alpha,FDRflag);
                    fprintf(fid,'%s,',str);
                    for j = 2:NCondSteps
                        if CondEffectsSIGN(j-1) ~= 0
                            fprintf(fid,'%7.3f*,',CondEffects(j-1));
                        else
                            fprintf(fid,'%7.3f,',CondEffects(j-1));
                        end
                    end
                    
                    fprintf(fid,'\n');
                end
                
            end
            
        case '4'
            if PrintHeaderFlag
                fprintf(fid,'Alpha = %s\n',StrAlpha);
                fprintf(fid,'\n%4s,%-30s,%-10s,%-10s,%7s,%7s,%7s,','ind','Region','Hemi','Measure','a','cP','b');
                fprintf(fid,'\n');
                PrintHeaderFlag = 0;
            end
            % Check to make sure both interactions are significant
            % Model 1
            a = getfield(P.Model1{1},[P.names.X]);
            % Model 2
            cP =getfield(P.Model2,[P.names.X]);
            b = getfield(P.Model2,[P.names.M{1} '1']);
            % Check the significance of the conditional effect
            Limits = getfield(P.AB1{1}.BCaci,AlphaLimit);
            if prod(Limits) > 0
                %fprintf(fid,'%7d,',1);
                %fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,num2str(i),'?','area');
                fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,outNames{i},Hemi{i},MeasurementType{i});
                str = subfnCreateParamStr(a,alpha,FDRflag);
                fprintf(fid,'%s,',str);
                str = subfnCreateParamStr(cP,alpha,FDRflag);
                fprintf(fid,'%s,',str);
                str = subfnCreateParamStr(b,alpha,FDRflag);
                fprintf(fid,'%s,',str);
                fprintf(fid,'\n');
                
            end
            % print out results
            
    end
    
end


function out = subfnCreateParamStr(input,alpha,FDRflag)
BigSteps = 10.^[5:-1:0];
SmallSteps = 1./10.^([1:5]);
BigSigDig = 0;
SmallSigDig = 0;
for i = 1:length(BigSteps)
    if round(input.se/BigSteps(i))
        BigSigDig = 7 - i;
        break
    end
end
for i = 1:length(SmallSteps)
    if round(input.se/SmallSteps(i))
        SmallSigDig = i;
        break
    end
end
Str = eval(sprintf('sprintf(''%%0.%df'',%s)',SmallSigDig,input.beta));
if (abs(input.p) < alpha) & (FDRflag)
    out = sprintf('%s**',Str);
elseif (abs(input.p) < 0.05) & (FDRflag)
    out = sprintf('%s*',Str);
elseif (abs(input.p) < alpha)
    out = sprintf('%s*',Str);
else
    out = sprintf('%s',Str);
end

% if input.p < alpha
%     if abs(input.beta) < 0.001
%         out = sprintf('%7.2e*',input.beta);
%     else
%         out = sprintf('%7.3f*',input.beta);
%     end
% else
%     if abs(input.beta) < 0.001
%         out = sprintf('%7.2e',input.beta);
%     else
%         out = sprintf('%7.3f',input.beta);
%     end
% end
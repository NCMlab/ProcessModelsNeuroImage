function subfnPrintFSResults_beta(SelectedPath,alpha,fid)
% write Fressurfer results to tables
switch nargin
    case 2
        fid = 1;
    case 1
        alpha = 0.05;
        fid = 1;
    case 0
        SelectedPath = spm_select(1,'dir');
        alpha = 0.05;
        fid = 1;
end

cd(SelectedPath)
if exist('AnalysisParameters.mat')
    load AnalysisParameters
    load Results_0001
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
StrAlpha = num2str(alpha);
Dot=findstr(StrAlpha,'.');
AlphaLimit = ['alpha' StrAlpha(Dot+1:end)];

fprintf(fid,'==========================================\n')
fprintf(fid,'%s\n',FileName);

FDRflag = 0;

PrintHeaderFlag = 1;
%% Find FDR rates
switch ModelNum
    case '58'
        % Get p values
        a = zeros(Nvoxels,1);
        p  = zeros(Nvoxels,1);
        w = zeros(Nvoxels,1);
        cP = zeros(Nvoxels,1);
        b = zeros(Nvoxels,1);
        q = zeros(Nvoxels,1);
        v = zeros(Nvoxels,1);
        for i = 1:Nvoxels
            P = Parameters{i};
            % Model 1
            Sa = getfield(P.Model1{1},[P.names.X]);
            a(i) = Sa.p;
            Sp = getfield(P.Model1{1},[P.names.W]);
            p(i) = Sp.p;
            Sw = getfield(P.Model1{1},[P.names.X '_x_' P.names.W]);
            w(i) = Sw.p;
            % Model 2
            ScP =getfield(P.Model2,[P.names.X]);
            cP(i) = ScP.p;
            Sb = getfield(P.Model2,[P.names.M{1}]);
            b(i) = Sb.p;
            Sq = getfield(P.Model2,[P.names.W]);
            q(i) = Sq.p;
            Sv = getfield(P.Model2,[P.names.M{1} '_x_' P.names.W]);
            v(i) = Sv.p;
        end
        FDRa = FDR(a,alpha);
        FDRp = FDR(p,alpha);
        FDRw = FDR(w,alpha);
        FDRcP = FDR(cP,alpha);
        FDRb = FDR(b,alpha);
        FDRq = FDR(q,alpha);
        FDRv = FDR(v,alpha);
        if isempty(FDRa)
            FDRa = 10^-10;
        end
        if isempty(FDRp)
            FDRp = 10^-10;
        end
        if isempty(FDRw)
            FDRw = 10^-10;
        end
        if isempty(FDRcP)
            FDRcP = 10^-10;
        end
        if isempty(FDRb)
            FDRb = 10^-10;
        end
        if isempty(FDRq)
            FDRq = 10^-10;
        end
        if isempty(FDRv)
            FDRv = 10^-10;
        end
        if FDRflag
            wTHR = FDRv;
        else
            wTHR = alpha;
        end
        if FDRflag
            vTHR = FDRv;
        else
            vTHR = alpha;
        end
    case '14'
        % Get p values
        a = zeros(Nvoxels,1);
        cP = zeros(Nvoxels,1);
        b = zeros(Nvoxels,1);
        q = zeros(Nvoxels,1);
        v = zeros(Nvoxels,1);
        for i = 1:Nvoxels
            P = Parameters{i};
            % Model 1
            Sa = getfield(P.Model1{1},[P.names.X]);
            a(i) = Sa.p;
            % Model 2
            ScP =getfield(P.Model2,[P.names.X]);
            cP(i) = ScP.p;
            Sb = getfield(P.Model2,[P.names.M{1}]);
            b(i) = Sb.p;
            Sq = getfield(P.Model2,[P.names.V]);
            q(i) = Sq.p;
            Sv = getfield(P.Model2,[P.names.M{1} '_x_' P.names.V]);
            v(i) = Sv.p;
        end
        FDRa = FDR(a,alpha);
        FDRcP = FDR(cP,alpha);
        FDRb = FDR(b,alpha);
        FDRq = FDR(q,alpha);
        FDRv = FDR(v,alpha);
        if isempty(FDRa)
            FDRa = 10^-10;
        end
        if isempty(FDRcP)
            FDRcP = 10^-10;
        end
        if isempty(FDRb)
            FDRb = 10^-10;
        end
        if isempty(FDRq)
            FDRq = 10^-10;
        end
        if isempty(FDRv)
            FDRv = 10^-10;
        end
        if FDRflag
            vTHR = FDRv;
        else
            vTHR = alpha;
        end
    case '7'
        % Get p values
        a = zeros(Nvoxels,1);
        p  = zeros(Nvoxels,1);
        w = zeros(Nvoxels,1);
        cP = zeros(Nvoxels,1);
        b = zeros(Nvoxels,1);
        for i = 1:Nvoxels
            P = Parameters{i};
            % Model 1
            Sa = getfield(P.Model1{1},[P.names.X]);
            a(i) = Sa.p;
            Sp = getfield(P.Model1{1},[P.names.W]);
            p(i) = Sp.p;
            Sw = getfield(P.Model1{1},[P.names.X '_x_' P.names.W]);
            w(i) = Sw.p;
            % Model 2
            ScP =getfield(P.Model2,[P.names.X]);
            cP(i) = ScP.p;
            Sb = getfield(P.Model2,[P.names.M{1}]);
            b(i) = Sb.p;
        end
        FDRa = FDR(a,alpha);
        FDRp = FDR(p,alpha);
        FDRw = FDR(w,alpha);
        FDRcP = FDR(cP,alpha);
        FDRb = FDR(b,alpha);
        if isempty(FDRa)
            FDRa = 10^-10;
        end
        if isempty(FDRp)
            FDRp = 10^-10;
        end
        if isempty(FDRw)
            FDRw = 10^-10;
        end
        if isempty(FDRcP)
            FDRcP = 10^-10;
        end
        if isempty(FDRb)
            FDRb = 10^-10;
        end
        if FDRflag
            wTHR = FDRw;
        else
            wTHR = alpha;
        end
    case '4'
        % Get p values
        a = zeros(Nvoxels,1);
        cP = zeros(Nvoxels,1);
        b = zeros(Nvoxels,1);
        for i = 1:Nvoxels
            P = Parameters{i};
            % Model 1
            Sa = getfield(P.Model1{1},[P.names.X]);
            a(i) = Sa.p;
            % Model 2
            ScP =getfield(P.Model2,[P.names.X]);
            cP(i) = ScP.p;
            Sb = getfield(P.Model2,[P.names.M{1} '1']);
            b(i) = Sb.p;
        end
        FDRa = FDR(a,alpha);
        FDRcP = FDR(cP,alpha);
        FDRb = FDR(b,alpha);
        if isempty(FDRa)
            FDRa = 10^-10;
        end
        if isempty(FDRcP)
            FDRcP = 10^-10;
        end
        if isempty(FDRb)
            FDRb = 10^-10;
        end
end



for i = 1:AnalysisParameters.Nvoxels
    P = Parameters{i};
    switch ModelNum
        case '58'
            NCondSteps = length(Parameters{1}.CondAB1);
            if PrintHeaderFlag
                fprintf(fid,'Interaction threshold = %0.4f, %0.4f\n',wTHR,vTHR);
                fprintf(fid,'%4s,%-30s,%-10s,%-10s,%7s,%7s,%7s,%7s,%7s,%7s,%7s,','ind','Region','Hemi','Measure','a','p','w','cP','b','q','v');
                for j = 2:NCondSteps
                    fprintf(fid,'%7.3f,',Parameters{1}.CondAB1{j}.probeValue);
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
                    fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,outNames{i},Hemi{i},MeasurementType{i});
                    str = subfnCreateParamStr(a,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(p,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(w,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(cP,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(b,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(q,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(v,alpha);
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
            NCondSteps = length(Parameters{1}.CondAB1);
            if PrintHeaderFlag
                fprintf(fid,'Interaction threshold = %0.4f\n',vTHR);
                fprintf(fid,'%4s,%-30s,%-10s,%-10s,%7s,%7s,%7s,%7s,%7s,','ind','Region','Hemi','Measure','a','cP','b','q','v');
                for j = 2:NCondSteps
                    fprintf(fid,'%7.3f,',Parameters{1}.CondAB1{j}.probeValue);
                end
                fprintf(fid,'\n');
                PrintHeaderFlag = 0;
            end
            % Check to make sure both interactions are significant
            % Model 1
            a = getfield(P.Model1{1},[P.names.X]);
            % Model 2
            cP =getfield(P.Model2,[P.names.X]);
            b = getfield(P.Model2,[P.names.M{1}]);
            q = getfield(P.Model2,[P.names.V]);
            v = getfield(P.Model2,[P.names.M{1} '_x_' P.names.V]);
            % Is the interaction significant?
            if abs(v.p) < vTHR
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
                    fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,outNames{i},Hemi{i},MeasurementType{i});
                    str = subfnCreateParamStr(a,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(cP,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(b,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(q,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(v,alpha);
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
        case '7'
            NCondSteps = length(Parameters{1}.CondAB1);
            if PrintHeaderFlag
                fprintf(fid,'Interaction threshold = %0.4f\n',wTHR);
                fprintf(fid,'%4s,%-30s,%-10s,%-10s,%7s,%7s,%7s,%7s,%7s,','ind','Region','Hemi','Measure','a','p','w','cP','b');
                for j = 2:NCondSteps
                    fprintf(fid,'%7.3f,',Parameters{1}.CondAB1{j}.probeValue);
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

                if ~isempty(find(CondEffects~=0))
                    % print out results
                    fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,outNames{i},Hemi{i},MeasurementType{i});
                    str = subfnCreateParamStr(a,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(p,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(w,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(cP,alpha);
                    fprintf(fid,'%s,',str);
                    str = subfnCreateParamStr(b,alpha);
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
                fprintf(fid,'%4d,%-30s,%-10s,%-10s,',i,outNames{i},Hemi{i},MeasurementType{i});
                str = subfnCreateParamStr(a,alpha);
                fprintf(fid,'%s,',str);
                str = subfnCreateParamStr(cP,alpha);
                fprintf(fid,'%s,',str);
                str = subfnCreateParamStr(b,alpha);
                fprintf(fid,'%s,',str);
                fprintf(fid,'\n');
                
            end
            % print out results
            
    end
    
end


function out = subfnCreateParamStr(input,alpha);
if input.p < alpha
    out = sprintf('%7.3f*',input.beta);
else
    out = sprintf('%7.3f',input.beta);
end
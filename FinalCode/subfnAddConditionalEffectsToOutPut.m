function OutData = subfnAddConditionalEffectsToOutPut(OutData,AnalysisParameters)
index = length(OutData) + 1;
ModelNum = AnalysisParameters.ModelNum;
Nmed = AnalysisParameters.Nmed;
Nthr = length(AnalysisParameters.Thresholds);
Thresholds = AnalysisParameters.Thresholds;
Nvoxels = AnalysisParameters.Nvoxels;
switch ModelNum
    case '4'
        % Conditional Effects
        for j = 1:Nmed
            for i = 1:Nthr
                thrStr = num2str(num2str(Thresholds(i)));
                OutData{index}.name = ['ABMed' num2str(j) 'sign_' num2str(Thresholds(i))];
                OutData{index}.data = zeros(Nvoxels,1);
                OutData{index}.field = ['AB' num2str(j) '{' num2str(j) '}.BCaci.alpha' thrStr(3:end)];
                OutData{index}.dataType = 2;
                index = index + 1;
            end
        end
        OutData{index}.name = 'k2';
        OutData{index}.data = zeros(Nvoxels,1);
        OutData{index}.field = ['k2.pointEst'];
        OutData{index}.dataType = 16;
        index = index + 1;
        
    case '7'
        for j = 1:Nmed
            for i = 1:Nthr
                thrStr = num2str(num2str(Thresholds(i)));
                OutData{index}.name = sprintf('CondABMed%d_pV000_sign%0.4f',j,Thresholds(i));
                OutData{index}.data = zeros(Nvoxels,1);
                OutData{index}.field = ['CondAB' num2str(j) '{1}.BCaci.alpha' thrStr(3:end)];
                OutData{index}.dataType = 2;
                index = index + 1;
            end
            OutData{index}.name = sprintf('CondABMed%d_pV000_pointEst',j);
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['CondAB' num2str(j) '{1}.pointEst'];
            OutData{index}.dataType = 16;
            index = index + 1;
        end
        
    case '14'
        for j = 1:Nmed
            for i = 1:Nthr
                thrStr = num2str(num2str(Thresholds(i)));
                OutData{index}.name = sprintf('CondABMed%d_pV000_sign%0.4f',j,Thresholds(i));
                OutData{index}.data = zeros(Nvoxels,1);
                OutData{index}.field = ['CondAB' num2str(j) '{1}.BCaci.alpha' thrStr(3:end)];
                OutData{index}.dataType = 2;
                index = index + 1;
            end
            OutData{index}.name = sprintf('CondABMed%d_pV000_pointEst',j);
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['CondAB' num2str(j) '{1}.pointEst'];
            OutData{index}.dataType = 16;
            index = index + 1;
        end
%         for j = 1:Nmed
%             OutData{index}.name = ['CondAB' num2str(j)];
%             OutData{index}.data = zeros(Nvoxels,1);
%             OutData{index}.field = ['CondAB' num2str(j) '{' num2str(j) '}.pointEst'];
%             index = index + 1;
% %             OutData{index}.name = ['CondABpValue' num2str(j)];
% %             OutData{index}.data = zeros(Nvoxels,1);
% %             OutData{index}.field = ['CondAB' num2str(j) '{' num2str(j) '}.pValue'];
% %             index = index + 1;
%             for i = 1:Nthr
%                 thrStr = num2str(num2str(Thresholds(i)));
%                 OutData{index}.name = ['CondABMed' num2str(j) 'sign_' num2str(Thresholds(i))];
%                 OutData{index}.data = zeros(Nvoxels,1);
%                 OutData{index}.field = ['CondAB' num2str(j) '{' num2str(j) '}.BCaci.alpha' thrStr(3:end)];
%                 index = index + 1;
%             end
%         end
        
    case '58'
        for j = 1:Nmed
            for i = 1:Nthr
                thrStr = num2str(num2str(Thresholds(i)));
                OutData{index}.name = sprintf('CondABMed%d_pV000_sign%0.4f',j,Thresholds(i));
                OutData{index}.data = zeros(Nvoxels,1);
                OutData{index}.field = ['CondAB' num2str(j) '{1}.BCaci.alpha' thrStr(3:end)];
                OutData{index}.dataType = 2;
                index = index + 1;
            end
            OutData{index}.name = sprintf('CondABMed%d_pV000_pointEst',j);
            OutData{index}.data = zeros(Nvoxels,1);
            OutData{index}.field = ['CondAB' num2str(j) '{1}.pointEst'];
            OutData{index}.dataType = 16;
            index = index + 1;
        end
        
    case '74'
        % Conditional Effects
        NProbe = length(AllParameters{1}.CondAB1);
        for j = 1:Nmed
            for k = 1:NProbe
                probeValue = AllParameters{1}.CondAB1{k}.probeValue;
                for i = 1:Nthr
                    thrStr = num2str(num2str(Thresholds(i)));
                    OutData{index}.name = sprintf('CondABMed%d_pV%0.2f_sign%0.4f',j,probeValue,Thresholds(i));
                    OutData{index}.data = zeros(Nvoxels,1);
                    OutData{index}.field = ['CondAB' num2str(j) '{' num2str(k) '}.BCaci.alpha' thrStr(3:end)];
                    OutData{index}.dataType = 2;
                    index = index + 1;
                    
                end
            end
        end
end
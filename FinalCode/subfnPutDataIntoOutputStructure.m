function OutData = subfnPutDataIntoOutputStructure(OutData,AllParameters,AnalysisParameters)
% put the data into the output structure
Nvoxels = AnalysisParameters.Nvoxels;
OutDataForProbesFlag = 0;
for i = 1:Nvoxels
    if mod(i,1000) == 0
        fprintf(1,'Working on voxel: %d of %d\n',i,Nvoxels);
    end

    if ~isempty(AllParameters{i})
        VoxelIndices(i) = 1;
        ImageVoxelIndices(i) = AnalysisParameters.Indices(i);
        % this limit for the length of OutData maye change due to the
        % discovery of significant probe values in the results. Therefore, %
        % it must remain as a check of the length.
        if ~OutDataForProbesFlag
            for j = 1:length(OutData)
                % If an OutData structure refers to something with a probe
                % value then there needs to be a loop here that cycles over all
                % probe values and writes out the images. This needs to be done
                % here because we cannot assume that each voxel will have the
                % same moderating variable, i.e. there may be a voxelwise
                % moderator.
                
                % only add the new structures for the probe values if they
                % have NOT been made before.
                if ~isempty(strfind(OutData{j}.name,'000'))
                    
                    % This is a probe measure.
                    % Since this is a probe measure, the OutData cell is
                    % replaced with a structure, replicating itself for each of
                    % the probeValues.
                    % How many probe values are there for this voxel?
                    Nprobe = length(AllParameters{i}.CondAB1);
                    if Nprobe > 1
                        OutDataForProbesFlag = 1;
                    end
                    probeValues = zeros(Nprobe,1);
                    for k = 1:Nprobe
                        probeValues(k) = AllParameters{i}.CondAB1{k}.probeValue;
                    end
                    tempOutData = {};
                    index = 1;
                    tempName = OutData{j}.name;
                    tempField = OutData{j}.field;
                    for k = 1:Nprobe
                        tempOutData{index}.name = strrep(tempName,'000',sprintf('%0.2f',probeValues(k)));
                        tempOutData{index}.field = strrep(tempField,'000',num2str(k));
                        tempOutData{index}.data = zeros(Nvoxels,1);
                        tempOutData{index}.dataType = 2;
                        index = index + 1;
                    end
                   OutData{j} = [];
                    index = length(OutData);
                    for k = 1:Nprobe
                        OutData{index+k} = tempOutData{k};
                    end
                end
            end
        end
         % remove the OutData cell and add the new structure to
                    % the end.
                    
                OutData = OutData(~cellfun('isempty',OutData));
                    
        % This double checking may be slow, but it is done so that if you
        % are in the first voxel with probe values that the new additions
        % ot the OutData are still cycled over.
        for j = 1:length(OutData)
            % the confidence intervals need to be treated special
            if strfind(OutData{j}.field,'alpha')
                % This try/catch block is to see if there is more then one
                % probe value. If not then there is an error and it is
                % skipped. But this assumes that we know that there are
                % saved probe values.
                try
                    thresholds = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                    if prod(thresholds)>0
                        OutData{j}.data(i) = 1;
                        OutData{j}.dataType = 2;
                    end
                catch
                end
                
            else
                % fprintf(1,'%s\n', ['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                % fprintf(1,'%s\n',['AllParameters{' num2str(i) '}.' OutData{j}.field]);
                OutData{j}.data(i) = eval(['AllParameters{' num2str(i) '}.' OutData{j}.field]);
            end
        end
    end
end

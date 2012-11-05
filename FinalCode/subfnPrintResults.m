function subfnPrintResults(Parameters,fid)
if nargin == 1
    fid = 1;
end
ModelNum = Parameters.ModelNum;
switch ModelNum
    case '1'
        fprintf(fid,'======================================================\n');
        fprintf(fid,'Model = %s\n',Parameters.ModelNum);
        fprintf(fid,'\tY = %s\n',Parameters.Yname);
        fprintf(fid,'\tX = %s\n',Parameters.Xname);
        fprintf(fid,'\tM = %s\n',Parameters.Mname);
        fprintf(fid,'Sample size = %d\n\n',Parameters.SampleSize);
        
        fprintf(fid,'******************************************************\n')
        fprintf(fid,'Outcome: %s\n\n',Parameters.Model1.Outcome)
        subfnPrintModelSummary(Parameters.Model1.Model,fid)
        subfnPrintModelResults(Parameters.Model1,fid)
        fprintf(fid,'\nR-square increase due to interaction:\n');
        subfnPrintModelSummary(Parameters.DiffModel,fid)
        fprintf(fid,'Moderator value(s) defining Johnson-Neyman significance region(s)\n');
        fprintf(fid,'\t%8.4f\n',Parameters.JNvalue)
        fprintf(fid,'Conditional effect of %s on %s at values of the moderator (%s):\n',...
            Parameters.Xname,Parameters.Yname,Parameters.Mname);
        fprintf(fid,'%8s\t%8s\t%8s\t%8s\t%8s\n',...
            Parameters.Mname,'Effect','boot se','lowerCI','upperCI')
        for k = 2:length(Parameters.CondMod)
            fprintf(fid,'%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n',...
                Parameters.CondMod{k}.probeValue,...
                Parameters.CondMod{k}.pointEst,...
                Parameters.CondMod{k}.bootSE,...
                Parameters.CondMod{k}.BCaci.alpha05(1),...
                Parameters.CondMod{k}.BCaci.alpha05(2));
        end
    case '4'
        % MODEL 4
        fprintf(1,'======================================================\n');
        fprintf(1,'Model = %s\n',Parameters.ModelNum);
        fprintf(1,'\tY = %s\n',Parameters.Yname);
        fprintf(1,'\tX = %s\n',Parameters.Xname);
        fprintf(1,'\tM = %s\n',Parameters.Mname);
        fprintf(1,'Sample size = %d\n\n',Parameters.SampleSize);
        
        
        fprintf(1,'Indirect effect of %s on %s via %s (a*b pathway)\n',Parameters.Xname,Parameters.Yname,Parameters.Mname);
        fprintf(1,'%8s\t%8s\t%8s\t%8s\n','Effect','Boot SE','BootLLCI','BootULCI');
        fprintf(1,'%8.4f\t%8.4f\t%8.4f\t%8.4f\n',Parameters.AB1{1}.pointEst,Parameters.AB1{1}.bootSE,Parameters.AB1{1}.BCaci.alpha05(1),Parameters.AB1{1}.BCaci.alpha05(2));
        fprintf(1,'Preacher and Kelley (2011) Kappa-squared\n');
        fprintf(1,'%8s\t%8s\t%8s\t%8s\n','Effect','Boot SE','BootLLCI','BootULCI');
        fprintf(1,'%8.4f\t%8.4f\t%8.4f\t%8.4f\n',Parameters.k2.pointEst,-99,Parameters.k2.PERci{1}.alpha05(1),Parameters.k2.PERci{1}.alpha05(2));
        % Print out Model 1
        for i = 1:size(Parameters.Model1,2)
            fprintf(fid,'******************************************************\n')
            fprintf(1,'Outcome: %s\n\n',Parameters.Model1{i}.Outcome)
            subfnPrintModelSummary(Parameters.Model1{i}.Model,fid)
            subfnPrintModelResults(Parameters.Model1{i},fid)
        end
        
        % Print out Model 2
        fprintf(fid,'******************************************************\n')
        fprintf(1,'Outcome: %s\n\n',Parameters.Model2.Outcome)
        subfnPrintModelSummary(Parameters.Model2.Model,fid)
        subfnPrintModelResults(Parameters.Model2,fid)
        
        % Print out Model 3
        fprintf(fid,'******************************************************\n')
        fprintf(1,'Outcome: %s\n\n',Parameters.Model3.Outcome)
        subfnPrintModelSummary(Parameters.Model3.Model,fid)
        subfnPrintModelResults(Parameters.Model3,fid)
    case '7'
        % MODEL 7
        %
        %  W  M
        %   \/ \
        %   /   \
        %  X     Y
        %
        fprintf(1,'======================================================\n');
        fprintf(1,'Model = %s\n',Parameters.ModelNum);
        fprintf(1,'\tY = %s\n',Parameters.Yname);
        fprintf(1,'\tX = %s\n',Parameters.Xname);
        fprintf(1,'\tM = %s\n',Parameters.Mname);
        fprintf(1,'\tW = %s\n',Parameters.Wname);
        fprintf(1,'Sample size = %d\n\n',Parameters.SampleSize);
        % Print out Model 1
        for i = 1:size(Parameters.Model1,2)
            fprintf(fid,'******************************************************\n')
            fprintf(1,'Outcome: %s\n\n',Parameters.Model1{i}.Outcome)
            subfnPrintModelSummary(Parameters.Model1{i}.Model,fid)
            subfnPrintModelResults(Parameters.Model1{i},fid)
        end
        % Print out Model 2
        fprintf(fid,'******************************************************\n')
        fprintf(1,'Outcome: %s\n\n',Parameters.Model2.Outcome)
        subfnPrintModelSummary(Parameters.Model2.Model,fid)
        subfnPrintModelResults(Parameters.Model2,fid)
        % Print out Model 3
        fprintf(fid,'******************************************************\n')
        fprintf(1,'Outcome: %s\n\n',Parameters.Model3.Outcome)
        subfnPrintModelSummary(Parameters.Model3.Model,fid)
        subfnPrintModelResults(Parameters.Model3,fid)
        for i = 1:length(Parameters.Model1)
            fprintf(fid,'Conditional effect of %s on %s at values of the moderator (%s):\n',...
                Parameters.Xname,Parameters.Yname,Parameters.Mname);
            fprintf(fid,'%8s\t%8s\t%8s\t%8s\t%8s\n',...
                Parameters.Mname,'Effect','boot se','lowerCI','upperCI')
                 ThresholdField = fields(Parameters.CondAB1{1}.BCaci);
                ThresholdField = ThresholdField{1};           
            for k = 2:length(Parameters.CondAB1)


                limits = getfield(Parameters.CondAB1{k}.BCaci,ThresholdField);
                fprintf(fid,'%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n',...
                    Parameters.CondAB1{k}.probeValue,...
                    Parameters.CondAB1{k}.pointEst,...
                    Parameters.CondAB1{k}.bootSE,...
                    limits(1),...
                    limits(2));
            end
        end
        
end
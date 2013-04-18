function selected_PCs_LOOCV = subfnCompileClusterLOOCV(JobsDir,NPCs,NSub)
% Findthe LOOCV results files
F = dir(fullfile(JobsDir,'LOOCV*.mat'));
% Create matrix of all possible cominations of PCs
combo_matrix = boolean_enumeration_f(NPCs);
% how many combos are there?
NCombos = size(combo_matrix,1);

LOOCVMatrix = zeros(NSub,NCombos);
for i = 1:NSub
    temp = load(fullfile(JobsDir,F(i).name));
    LOOCVMatrix(i,:) = temp.LOOmediationMatrix;
end
clear temp
sLOOCVMatrix = sum(abs(LOOCVMatrix),1);
selected_PCs_LOOCV = find(combo_matrix(find(sLOOCVMatrix == max(sLOOCVMatrix)),:));
% Clean up the files created
for i = 1:NSub 
    unix(sprintf('rm %s',fullfile(JobsDir,F(i).name)));
end
unix(['rm job_LOOCV.sh*']);
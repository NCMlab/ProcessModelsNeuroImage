function LOOerrorMatrix = subfnCalcLOOAllModels(LOOdata,Totaldata,NPCs,SubjectIndex)
remove_row_means = 1;
combo_matrix = boolean_enumeration_f(NPCs);
NCombos = size(combo_matrix,1);

% apply PCA
    % but squeeze out the multiple mediator dimension
    [lambdas, eigenimages_noZeroes, w] = pca_f(LOOdata.M', remove_row_means);
    tempssf = squeeze(LOOdata.M) * eigenimages_noZeroes;
    tempssfSubset = tempssf(:,1:NPCs);
    % Fit the regression model with all PCs
    tempdata = LOOdata;
    tempdata.M = tempssfSubset;
  

    for j = 1:NCombos
        selected_PCs = find(combo_matrix(j,:));
        tempdata = LOOdata;
        tempdata.M = tempssfSubset(:,selected_PCs);
        %behav_fit_coef = FullModelbehav_fit_coef([1 selected_PCs+1 NPCs+2:end]);
        
        behav_fit_coef = subfnCallRegressPCs(tempdata,LOOdata.ModelNum);
        

        % create the SSF image
        temp = eigenimages_noZeroes(:, selected_PCs) * behav_fit_coef(1 + 1:1 + length(selected_PCs));  %nuisance regressors stay silent
        behav_fit_composite_PC_image = temp / norm(temp);
        %%%%% Obtain SSFs associated with the normalized best linear behavioral-fit image %%%%%
        behav_fit_composite_PC_image_ssfs = squeeze(LOOdata.M) * behav_fit_composite_PC_image;
        % refit the regression model 
        behav_fit_coef2 = subfnregress(LOOdata.Y,[behav_fit_composite_PC_image_ssfs LOOdata.X]);
        % forward apply this SSF to the left out subjects raw data
        % predict the left out subject
        predictedM = squeeze(Totaldata.M(SubjectIndex,:,:))'*behav_fit_composite_PC_image;
        
        % predict the left out subject
        predictedY = subfnLOOPredictPCs(Totaldata,behav_fit_coef2,LOOdata.ModelNum,predictedM,length(selected_PCs),SubjectIndex);
        LOOerrorMatrix(j) = (predictedY - Totaldata.Y(SubjectIndex))^2;
    end
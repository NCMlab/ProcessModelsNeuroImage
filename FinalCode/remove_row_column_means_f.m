function matrix=remove_row_column_means_f(data);
%function matrix=remove_row_column_means_f(data);
%
%
%
% Zero means DATA row and column wise

dim=size(data);
matrix=data-ones(dim(1),1)*mean(data,1)-mean(data,2)*ones(1,dim(2))+ones(dim(1),dim(2))*mean(mean(data));


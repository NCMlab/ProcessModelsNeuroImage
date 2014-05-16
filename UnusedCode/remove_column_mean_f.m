function matrix=remove_column_mean_f(data);
%function matrix=remove_column_mean_f(data);
%
%
%
% Zero means column wise

dim=size(data);
matrix=data-ones(dim(1),1)*mean(data,1);

